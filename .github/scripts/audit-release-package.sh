#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 PACKAGE_DIRECTORY" >&2
  exit 2
fi

package_dir=$(readlink -f "$1")

fail() {
  echo "release audit: $*" >&2
  exit 1
}

require_file() {
  [[ -f "$package_dir/$1" ]] || fail "missing file: $1"
}

require_executable() {
  [[ -x "$package_dir/$1" ]] || fail "not executable: $1"
}

version_at_most() {
  local actual=$1
  local maximum=$2
  [[ -z "$actual" ]] ||
    [[ $(printf '%s\n%s\n' "$actual" "$maximum" | sort -Vu | tail -n 1) == "$maximum" ]]
}

for command in file ldd objdump readelf readlink sort strings; do
  command -v "$command" >/dev/null || fail "required audit tool not found: $command"
done

for path in \
  README.md BENCHMARKS.md BUILD_INFO.txt LICENSE NOTICE \
  bin/megahit bin/megahit_core bin/megahit_core_popcnt \
  bin/megahit_core_no_hw_accel \
  lib/libdeflate.so.0 lib/libgcc_s.so.1 lib/libstdc++.so.6 lib/libz.so.1
do
  require_file "$path"
done

for path in \
  megahit bin/megahit bin/megahit_core bin/megahit_core_popcnt \
  bin/megahit_core_no_hw_accel bin/megahit_toolkit
do
  require_executable "$path"
done

[[ $(readlink "$package_dir/megahit") == bin/megahit ]] ||
  fail "top-level megahit must be a relative link to bin/megahit"
[[ $(readlink "$package_dir/bin/megahit_toolkit") == megahit_core_no_hw_accel ]] ||
  fail "megahit_toolkit must be a relative link to the portable core"
[[ ! -e "$package_dir/bin/rabbitma" ]] ||
  fail "unexpected second public CLI: bin/rabbitma"
[[ $(sed -n '1p' "$package_dir/bin/megahit") == '#!/usr/bin/env python3' ]] ||
  fail "megahit must use the portable python3 env shebang"

while IFS= read -r -d '' link; do
  resolved=$(readlink -f "$link")
  case "$resolved" in
    "$package_dir"/*) ;;
    *) fail "symlink escapes the package: ${link#"$package_dir"/}" ;;
  esac
done < <(find "$package_dir" -type l -print0)

bundled_libraries='libz.so.1 libdeflate.so.0 libstdc++.so.6 libgcc_s.so.1'

while IFS= read -r -d '' elf; do
  readelf -h "$elf" >/dev/null 2>&1 || continue
  relative=${elf#"$package_dir"/}

  machine=$(readelf -h "$elf" | awk -F: '/Machine:/ { sub(/^[[:space:]]+/, "", $2); print $2; exit }')
  [[ "$machine" == *X86-64* ]] || fail "$relative has unexpected machine type: $machine"

  while IFS= read -r needed; do
    case "$needed" in
      libc.so.6|libdl.so.2|libm.so.6|libpthread.so.0|ld-linux-x86-64.so.2) ;;
      libz.so.1|libdeflate.so.0|libstdc++.so.6|libgcc_s.so.1)
        [[ -f "$package_dir/lib/$needed" ]] ||
          fail "$relative needs an unbundled library: $needed"
        ;;
      *) fail "$relative has an unapproved DT_NEEDED dependency: $needed" ;;
    esac
  done < <(readelf -d "$elf" | sed -n 's/^.*Shared library: \[\(.*\)\]$/\1/p')

  max_glibc=$(readelf --version-info "$elf" 2>/dev/null |
    grep -oE 'GLIBC_[0-9.]+' | sort -Vu | tail -n 1 || true)
  version_at_most "$max_glibc" GLIBC_2.17 ||
    fail "$relative requires $max_glibc (maximum is GLIBC_2.17)"

  max_glibcxx=$(readelf --version-info "$elf" 2>/dev/null |
    grep -oE 'GLIBCXX_[0-9.]+' | sort -Vu | tail -n 1 || true)
  version_at_most "$max_glibcxx" GLIBCXX_3.4.19 ||
    fail "$relative requires $max_glibcxx (maximum is GLIBCXX_3.4.19)"

  max_cxxabi=$(readelf --version-info "$elf" 2>/dev/null |
    grep -oE 'CXXABI_[0-9.]+' | sort -Vu | tail -n 1 || true)
  version_at_most "$max_cxxabi" CXXABI_1.3.7 ||
    fail "$relative requires $max_cxxabi (maximum is CXXABI_1.3.7)"

  if strings "$elf" | grep -Eq '/home/|/usr/local|/opt/rabbitma-deps'; then
    fail "machine-local build path embedded in $relative"
  fi

  printf '%s: GLIBC=%s GLIBCXX=%s CXXABI=%s\n' \
    "$relative" "${max_glibc:-none}" "${max_glibcxx:-none}" "${max_cxxabi:-none}"
done < <(find "$package_dir/bin" "$package_dir/lib" -type f -print0)

for core in megahit_core megahit_core_popcnt megahit_core_no_hw_accel; do
  elf="$package_dir/bin/$core"
  interpreter=$(readelf -l "$elf" |
    sed -n 's/^.*Requesting program interpreter: \(.*\)]$/\1/p')
  [[ "$interpreter" == /lib64/ld-linux-x86-64.so.2 ]] ||
    fail "$core has unexpected interpreter: $interpreter"
  readelf -d "$elf" |
    grep -Eq '\((RPATH|RUNPATH)\).*[[]\$ORIGIN/\.\./lib[]]' ||
    fail "$core does not search its package-local lib directory"
  if env -u LD_LIBRARY_PATH ldd "$elf" | grep -q 'not found'; then
    fail "$core has an unresolved dynamic dependency"
  fi
done

# REP-prefixed BSF can be displayed as TZCNT while remaining compatible on
# pre-BMI CPUs when its operand is nonzero.  The instructions below do not
# have that compatibility property and must never enter the portable core.
if objdump -d "$package_dir/bin/megahit_core_no_hw_accel" |
    grep -Eq '\b(popcnt|pext|pdep|bzhi|shlx|shrx|sarx|andn|blsi|blsr|blsmsk|lzcnt|mulx)\b'
then
  fail "portable core contains an instruction requiring POPCNT/BMI"
fi

for library in $bundled_libraries; do
  [[ -f "$package_dir/lib/$library" ]] || fail "missing bundled $library"
done

echo "release audit: package structure, ABI, ISA and dependency closure passed"
