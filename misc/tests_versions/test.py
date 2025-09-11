#!/usr/bin/env python3
import subprocess, os, time
from pathlib import Path
import argparse
import sys
import csv
import json

STEPS = ["configure", "make", "examples", "tests"]


def get_arguments():
    parser = argparse.ArgumentParser(description="Test multiple OCaml and Sundials versions for the OCaml port.")
    parser.add_argument("--ocaml", "-o", nargs="+", help="List of OCaml switch versions to test", dest="ocaml_versions")
    parser.add_argument("--sundials", "-s", nargs="+", help="List of Sundials versions to test", dest="sundials_versions")
    parser.add_argument("--flag", "-f", default="", help="Optional flag to pass to the ./configure script")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")

    log_group = parser.add_mutually_exclusive_group()
    log_group.add_argument("--logs", "-l", metavar="LOGPATH", help="Create a directory for all logs")
    log_group.add_argument("--shortlogs", metavar="LOGPATH", help="Use single-file (short) logging")

    parser.add_argument("--silent", action="store_true", help="Silent mode: only show a progress bar")
    parser.add_argument("--step", choices=STEPS, default="make", help="Stop after this step (default: tests)")
    parser.add_argument("--sundials-dir", "-d", default="sundials", help="Directory containing Sundials versions")
    parser.add_argument("--sundialsml-dir", "-m", default="sundialsml", help="Directory containing the OCaml port")
    return parser.parse_args()


def print_progress_bar(current, total, stage, ocaml, sundials, length=35, _max=[0]):
    progress = int(length * current // total)
    bar = '[' + '#' * progress + '-' * (length - progress) + ']'
    percent = int(100 * current / total)
    message = f'{bar} {percent:3d}% : {stage} ({ocaml}/{sundials})'
    if len(message) > _max[0]:
        _max[0] = len(message)
    message += ' ' * (_max[0] - len(message))
    sys.stdout.write('\r' + message)
    sys.stdout.flush()
    if current == total:
        print()
        _max[0] = 0


def clean_dir(base_dir, sundialsml_dir, silent=False, log_file=None):
    sundialsml_path = base_dir / sundialsml_dir
    cmd = ["make", "distclean"]

    cmd_str = f"Running: {' '.join(cmd)} (cwd: {sundialsml_path})"
    if log_file:
        log_file.write(f"{cmd_str}\n")
        log_file.flush()
    if not silent:
        print(cmd_str)

    kwargs = {}
    if silent and not log_file:
        kwargs['stdout'] = subprocess.DEVNULL
        kwargs['stderr'] = subprocess.DEVNULL
    elif log_file:
        kwargs['stdout'] = log_file
        kwargs['stderr'] = log_file
    subprocess.run(cmd, cwd=sundialsml_path, check=False, **kwargs)


def run_make(cmd, cwd, silent, log_file):
    cmd_str = f"Running: {' '.join(cmd)} (cwd: {cwd})"
    if log_file:
        log_file.write(f"{cmd_str}\n")
        log_file.flush()
    if not silent:
        print(cmd_str)

    kwargs = {}
    if silent and not log_file:
        kwargs['stdout'] = subprocess.DEVNULL
        kwargs['stderr'] = subprocess.DEVNULL
    elif log_file:
        kwargs['stdout'] = log_file
        kwargs['stderr'] = log_file
    return subprocess.run(cmd, cwd=cwd, **kwargs)


def main():
    args = get_arguments()
    step_index = STEPS.index(args.step)

    log_file = open(args.shortlogs, "a") if args.shortlogs else None
    log_path = Path(args.logs) if args.logs else None

    if args.ocaml_versions:
        ocaml_versions = args.ocaml_versions
    else:
        cmd = ["opam", "switch", "list", "--short"]
        cmd_str = f"Running: {' '.join(cmd)}"
        if log_file:
            log_file.write(f"{cmd_str}\n")
            log_file.flush()
        if not args.silent:
            print(cmd_str)
        result = subprocess.run(cmd, capture_output=True, text=True)
        ocaml_versions = result.stdout.strip().splitlines()

    sundials_dir = Path(args.sundials_dir)
    if args.sundials_versions:
        sundials_versions = args.sundials_versions
    else:
        sundials_versions = sorted([
            entry for entry in os.listdir(sundials_dir)
            if (sundials_dir / entry / "install").is_dir()
        ])

    results = {ocaml: {} for ocaml in ocaml_versions}
    base_dir = Path.cwd()
    total = len(ocaml_versions) * len(sundials_versions)
    current = 0
    global_status = {
        "ocaml_version": ocaml_versions,
        "sundials_version": sundials_versions,
        "tests": [],
        "steps": STEPS
    }

    if log_path:
        log_path.mkdir(parents=True, exist_ok=True)
        status_json_path = log_path / "status.json"
        with open(status_json_path, "w") as status_file:
            json.dump(global_status, status_file, indent=2)

    for ocaml_ver in ocaml_versions:
        with open("result.csv", "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["OCaml\\Sundials"] + sundials_versions)
            for ocaml in ocaml_versions:
                row = [ocaml] + [results[ocaml].get(s, "") for s in sundials_versions]
                writer.writerow(row)
        for sundials_ver in sundials_versions:
            current += 1
            results[ocaml_ver][sundials_ver] = "success"
            key = f"{ocaml_ver}-{sundials_ver}" if log_path else None
            key_path = key_status_path = None
            key_status = {
                "state" : "running",
                "step" : {i: "waiting" for i in STEPS},
            }

            if log_path:
                key_path = log_path / key
                key_status_path = key_path / "status.json"
                key_path.mkdir(parents=True, exist_ok=True)
                with open(key_status_path, "w") as key_status_file:
                    json.dump(key_status, key_status_file, indent=2)

            job_log_path = None
            if log_path:
                job_log_path = log_path / f"{ocaml_ver}_{sundials_ver}"
                job_log_path.mkdir(parents=True, exist_ok=True)

            if log_file:
                log_file.write(f"\n--- {ocaml_ver}/{sundials_ver} ---\n")

            sundials_install_path = base_dir / args.sundials_dir / sundials_ver / "install"
            config_cmd = ["opam", "exec", "--switch", ocaml_ver, "--", "./configure", f"SUNDIALS_DIR={sundials_install_path}", f"CFLAGS=\"-Wno-error=dangling-pointer\""]

            if args.flag:
                config_cmd.append(args.flag)

            if args.silent:
                print_progress_bar(current, total, "configure", ocaml_ver, sundials_ver)

            if job_log_path:
                job_log_file_path = job_log_path / "configure.log"
                log_file = open(job_log_file_path, "w")

            make_proc = run_make(
                config_cmd, 
                base_dir / args.sundialsml_dir, 
                args.silent, 
                log_file
            )

            if job_log_path:
                log_file.close()

            if make_proc.returncode != 0:
                results[ocaml_ver][sundials_ver] = "configure"
                continue
            elif step_index == 0:
                continue

            if job_log_path:
                job_log_file_path = job_log_path / "make.log"
                log_file = open(job_log_file_path, "w")

            if args.silent:
                print_progress_bar(current, total, "make", ocaml_ver, sundials_ver)

            make_proc = run_make(
                ["opam", "exec", "--switch", ocaml_ver, "--", "make", "all", "-j"],
                base_dir / args.sundialsml_dir, args.silent, log_file)

            if job_log_path:
                log_file.close()

            if make_proc.returncode != 0:
                results[ocaml_ver][sundials_ver] = "make"
                continue
            elif step_index == 1:
                continue

            if job_log_path:
                job_log_file_path = job_log_path / "make_examples.log"
                log_file = open(job_log_file_path, "w")

            if args.silent:
                print_progress_bar(current, total, "make examples", ocaml_ver, sundials_ver)

            examples_path = base_dir / args.sundialsml_dir / "examples"
            make_ex = run_make(["opam", "exec", "--switch", ocaml_ver, "--", "make", "all", "-j"], examples_path, args.silent, log_file)

            if job_log_path:
                log_file.close()

            if make_ex.returncode != 0:
                results[ocaml_ver][sundials_ver] = "make examples"
                continue
            elif step_index == 2:
                continue

            if job_log_path:
                job_log_file_path = job_log_path / "exec_tests.log"
                log_file = open(job_log_file_path, "w")

            if args.silent:
                print_progress_bar(current, total, "tests", ocaml_ver, sundials_ver)
            make_test = run_make(["opam", "exec", "--switch", ocaml_ver, "--", "make", "tests.opt.log", "-j"], examples_path, args.silent, log_file)

            if make_test.returncode != 0:
                results[ocaml_ver][sundials_ver] = "exec tests"
            clean_dir(base_dir, args.sundialsml_dir, silent=True, log_file=log_file)

    if log_file:
        log_file.close()

    with open("result.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["OCaml\\Sundials"] + sundials_versions)
        for ocaml in ocaml_versions:
            row = [ocaml] + [results[ocaml].get(s, "") for s in sundials_versions]
            writer.writerow(row)

    print("\nRésultats sauvegardés dans result.csv")

if __name__ == "__main__":
    exit(main())
