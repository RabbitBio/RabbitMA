#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 PACKAGE_DIRECTORY" >&2
  exit 2
fi

package_dir=$(readlink -f "$1")
runtime_path=${RABBITMA_RUNTIME_PATH:-/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin}
work_dir=$(mktemp -d '/tmp/RabbitMA runtime.XXXXXX')
trap 'chmod -R u+w "$work_dir" 2>/dev/null || true; rm -rf "$work_dir"' EXIT

for command in bzip2 gzip ldd python3; do
  PATH=$runtime_path command -v "$command" >/dev/null || {
    echo "runtime test: missing command: $command" >&2
    exit 1
  }
done

PATH=$runtime_path python3 - \
  "$package_dir/bin/megahit" "$work_dir/megahit.pyc" <<'PY'
import py_compile
import sys

if sys.version_info < (3, 5):
    raise SystemExit("RabbitMA requires Python 3.5 or newer")
py_compile.compile(sys.argv[1], cfile=sys.argv[2], doraise=True)
print("python runtime: %s" % sys.version.split()[0])
PY

if { printf 'write probe\n' >"$package_dir/.rabbitma-write-probe"; } 2>/dev/null; then
  rm -f "$package_dir/.rabbitma-write-probe"
  echo "runtime test: package mount is writable; expected a read-only install" >&2
  exit 1
fi

run_clean() {
  env -i \
    PATH="$runtime_path" \
    HOME="$work_dir" \
    TMPDIR="$work_dir" \
    "$@"
}

for core in megahit_core megahit_core_popcnt megahit_core_no_hw_accel; do
  if env -u LD_LIBRARY_PATH ldd "$package_dir/bin/$core" | grep -q 'not found'; then
    echo "runtime test: unresolved dependency in $core" >&2
    exit 1
  fi
done

run_clean "$package_dir/megahit" --version | tee "$work_dir/version.log"
grep -Eq '^RabbitMA [0-9]+\.[0-9]+\.[0-9]+ \(MEGAHIT core v1\.2\.9\)$' \
  "$work_dir/version.log"

run_clean "$package_dir/megahit" --test -t 2 2>&1 |
  tee "$work_dir/automatic.log"
run_clean "$package_dir/megahit" --test --no-hw-accel -t 2 2>&1 |
  tee "$work_dir/portable.log"

grep -q 'ALL DONE' "$work_dir/automatic.log"
grep -q 'ALL DONE' "$work_dir/portable.log"
grep -q 'without POPCNT and BMI2 support' "$work_dir/portable.log"

automatic_stats=$(grep -E '[0-9]+ contigs, total [0-9]+ bp' \
  "$work_dir/automatic.log" | tail -n 1 | sed 's/^.* - //')
portable_stats=$(grep -E '[0-9]+ contigs, total [0-9]+ bp' \
  "$work_dir/portable.log" | tail -n 1 | sed 's/^.* - //')
[[ -n "$automatic_stats" && "$automatic_stats" == "$portable_stats" ]] || {
  echo "runtime test: automatic and portable core statistics differ" >&2
  echo "automatic: $automatic_stats" >&2
  echo "portable:  $portable_stats" >&2
  exit 1
}

echo "runtime test: read-only, clean-environment and portable-core tests passed"
