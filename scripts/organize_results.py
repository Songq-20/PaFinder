#!/usr/bin/env python3
"""Consolidate per-MAG result directories without altering analysis evidence.

Usage:
    python scripts/organize_results.py /path/to/results            # preview
    python scripts/organize_results.py /path/to/results --apply    # move

Never run while a PaFinder job is still writing to the output directory.
"""

import argparse
import csv
import shutil
from pathlib import Path


def organize(root: Path, apply: bool = False) -> tuple[int, int]:
    root = root.expanduser().resolve()
    summary = root / "summary.tsv"
    if not summary.is_file():
        raise FileNotFoundError(f"Expected run summary: {summary}")

    with summary.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or "mag_id" not in rows[0]:
        raise ValueError("summary.tsv must contain at least one row with a mag_id column")

    ids = list(dict.fromkeys(row["mag_id"].strip() for row in rows))
    if any(not mid or mid in (".", "..") or Path(mid).name != mid for mid in ids):
        raise ValueError("Unsafe or empty MAG ID in summary.tsv")

    moves = []
    already = 0
    for mag_id in ids:
        source = root / mag_id
        dest = root / "mags" / mag_id
        if source.is_dir() and dest.exists():
            raise FileExistsError(f"Both source and destination exist: {source}, {dest}")
        if source.is_dir():
            moves.append((source, dest))
        elif dest.is_dir():
            already += 1
        else:
            raise FileNotFoundError(f"No result directory for {mag_id}: {source} or {dest}")

    for source, dest in moves:
        print(f"{'MOVE' if apply else 'PREVIEW'}: {source} -> {dest}")
    if apply and moves:
        (root / "mags").mkdir(exist_ok=True)
        for source, dest in moves:
            shutil.move(str(source), str(dest))

    print(f"MAGs in summary: {len(ids)}; to move: {len(moves)}; already organized: {already}")
    if not apply and moves:
        print("No files were moved. Add --apply to organize the completed run.")
    return len(moves), already


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", type=Path, help="Completed PaFinder output directory")
    parser.add_argument("--apply", action="store_true", help="Move per-MAG folders into results_dir/mags/")
    args = parser.parse_args()
    organize(args.results_dir, args.apply)


if __name__ == "__main__":
    main()
