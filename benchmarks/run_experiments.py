#!/usr/bin/env python3
"""Run the two experiments sequentially, so their CPU and I/O loads do not overlap."""

from pathlib import Path
import subprocess
import sys

from common import parse_arguments


def run_phase(command):
    # Forward an interruption once; the phase owns and cleans up its job groups.
    process = subprocess.Popen(command, start_new_session=True)
    try:
        status = process.wait()
        if status:
            raise subprocess.CalledProcessError(status, command)
    finally:
        if process.poll() is None:
            try:
                process.terminate()
            except ProcessLookupError:
                pass
            process.wait()


def main():
    args = parse_arguments(__doc__)
    root = Path(__file__).resolve().parent
    commands = []
    for script in ('compare_release_24.py', 'compare_concurrency.py'):
        command = [sys.executable, str(root / script), '--config', args.config]
        # Validate both plans before committing resources to either experiment.
        preview = subprocess.run(command + ['--dry-run'], check=True,
                                 stdout=subprocess.PIPE, universal_newlines=True)
        if args.dry_run:
            print(preview.stdout, end='')
        commands.append(command)
    if not args.dry_run:
        for command in commands:
            run_phase(command)


if __name__ == '__main__':
    main()
