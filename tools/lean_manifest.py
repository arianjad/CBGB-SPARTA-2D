#!/usr/bin/env python3
"""Create and finish the provenance receipt owned by run_lean.sh."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_output(root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), *args], capture_output=True, text=True, timeout=30
    )
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "git command failed")
    return result.stdout.strip()


def mpi_identity() -> str | None:
    try:
        result = subprocess.run(
            ["mpirun", "--version"], capture_output=True, text=True, timeout=30
        )
    except OSError:
        return None
    lines = (result.stdout or result.stderr).splitlines()
    return lines[0].strip() if lines else None


def write_atomic(path: Path, value: dict) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def start(args: argparse.Namespace) -> None:
    run_dir = Path(args.run_dir).resolve()
    repo = Path(args.repo_root).resolve()
    commit = git_output(repo, "rev-parse", "HEAD")
    if len(commit) != 40:
        raise RuntimeError("repository commit is not a full Git object ID")
    inputs = []
    for record in args.input_record:
        source_name, staged_name, role = record.split("\t", 2)
        source = Path(source_name).resolve()
        staged = run_dir / staged_name
        inputs.append(
            {
                "source_path": source.relative_to(repo).as_posix(),
                "staged_path": staged_name,
                "role": role,
                "source_sha256": sha256(source),
                "staged_sha256": sha256(staged),
            }
        )
    solver_hash = sha256(Path(args.solver_exe))
    verification = "declared"
    receipt_path = None
    if args.solver_receipt:
        receipt_file = Path(args.solver_receipt).resolve()
        receipt = json.loads(receipt_file.read_text(encoding="utf-8"))
        if (receipt.get("executable_sha256") != solver_hash or
                receipt.get("commit") != args.solver_commit or
                receipt.get("tag") != args.solver_tag):
            raise RuntimeError("solver build receipt does not match executable or declared source")
        verification = "build-receipt"
        receipt_path = str(receipt_file)
    manifest = {
        "schema_version": 1,
        "run_id": args.run_id,
        "case_id": args.case_id,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "running",
        "command": args.command_arg,
        "working_directory": str(run_dir),
        "repository": {
            "root": str(repo),
            "commit": commit,
            "dirty": bool(git_output(repo, "status", "--porcelain")),
        },
        "solver": {
            "selector": args.solver_selector,
            "tag": args.solver_tag,
            "commit": args.solver_commit,
            "executable": str(Path(args.solver_exe).resolve()),
            "executable_sha256": solver_hash,
            "source_identity": verification,
            "build_receipt": receipt_path,
        },
        "platform": {"os": platform.platform(), "mpi": mpi_identity()},
        "runtime_controls": {
            "mpi_ranks": args.mpi_ranks,
            **dict(item.split("=", 1) for item in args.runtime_control),
        },
        "inputs": inputs,
        "outputs": [],
        "started_utc": args.started,
        "finished_utc": None,
        "failure": None,
    }
    write_atomic(run_dir / "manifest.json", manifest)


def finish(args: argparse.Namespace) -> None:
    run_dir = Path(args.run_dir).resolve()
    path = run_dir / "manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    staged_inputs = {item["staged_path"] for item in manifest["inputs"]}
    roles = {
        "run.log": "run-log",
        "rc.sentinel": "exit-status",
        "start.stamp": "start-time",
        "manifest.txt": "legacy-receipt",
    }
    outputs = []
    for item in sorted(p for p in run_dir.rglob("*") if p.is_file() and p != path):
        rel = item.relative_to(run_dir).as_posix()
        if rel in staged_inputs:
            continue
        outputs.append(
            {
                "path": rel,
                "role": roles.get(rel, "generated-output"),
                "bytes": item.stat().st_size,
                "sha256": sha256(item),
            }
        )
    manifest["outputs"] = outputs
    manifest["finished_utc"] = args.finished
    manifest["status"] = "complete" if args.exit_code == 0 else "failed"
    manifest["failure"] = (
        None
        if args.exit_code == 0
        else {"exit_code": args.exit_code, "summary": f"SPARTA exited with status {args.exit_code}"}
    )
    write_atomic(path, manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    begin = subparsers.add_parser("start")
    begin.add_argument("--run-dir", required=True)
    begin.add_argument("--run-id", required=True)
    begin.add_argument("--case-id", required=True)
    begin.add_argument("--repo-root", required=True)
    begin.add_argument("--solver-selector", required=True)
    begin.add_argument("--solver-tag", required=True)
    begin.add_argument("--solver-commit", required=True)
    begin.add_argument("--solver-exe", required=True)
    begin.add_argument("--solver-receipt")
    begin.add_argument("--mpi-ranks", type=int, required=True)
    begin.add_argument("--runtime-control", action="append", default=[])
    begin.add_argument("--started", required=True)
    begin.add_argument("--command-arg", action="append", default=[])
    begin.add_argument("--input-record", action="append", default=[])
    end = subparsers.add_parser("finish")
    end.add_argument("--run-dir", required=True)
    end.add_argument("--exit-code", type=int, required=True)
    end.add_argument("--finished", required=True)
    args = parser.parse_args()
    (start if args.mode == "start" else finish)(args)


if __name__ == "__main__":
    main()
