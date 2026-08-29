#!/usr/bin/env python3
# Author: Rensc
# Date: 2026-08-29
# Version: dev001
# Function: Validate and generate BigWig files through pyBigWig for bwTools compatibility checks
# Input: command, BigWig/bedGraph/chrom-size paths
# Output: process exit status and optional BigWig output

import math
import sys
from pathlib import Path


def fail(message):
    print(f"[pyBigWig compatibility] ERROR: {message}", file=sys.stderr)
    return 1


def read_bedgraph(path):
    rows = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"bedGraph line {line_no} has fewer than four fields")
            chrom = fields[0]
            start0 = int(fields[1])
            end0 = int(fields[2])
            value = float(fields[3])
            if start0 < 0 or end0 <= start0 or not math.isfinite(value):
                raise ValueError(f"bedGraph line {line_no} contains invalid coordinates or value")
            rows.append((chrom, start0, end0, value))
    return rows


def read_chrom_sizes(path):
    rows = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"chrom sizes line {line_no} has fewer than two fields")
            chrom = fields[0]
            length = int(fields[1])
            if length < 1:
                raise ValueError(f"chrom sizes line {line_no} contains invalid length")
            rows.append((chrom, length))
    return rows


def import_pybigwig():
    try:
        import pyBigWig  # type: ignore
    except Exception as exc:
        raise RuntimeError(f"pyBigWig is not importable: {exc}") from exc
    return pyBigWig


def validate_file(bigwig_path, bedgraph_path, chrom_sizes_path, tolerance=1e-5):
    pyBigWig = import_pybigwig()
    expected = read_bedgraph(bedgraph_path)
    chrom_sizes = dict(read_chrom_sizes(chrom_sizes_path))
    bw = pyBigWig.open(str(bigwig_path))
    if bw is None or not bw.isBigWig():
        return fail("input is not recognized as BigWig")
    try:
        observed_chroms = dict(bw.chroms())
        if observed_chroms != chrom_sizes:
            return fail(f"chromosome sizes differ: observed={observed_chroms!r} expected={chrom_sizes!r}")

        observed = []
        for chrom, length in chrom_sizes.items():
            intervals = bw.intervals(chrom, 0, length)
            if intervals is None:
                continue
            for start0, end0, value in intervals:
                observed.append((chrom, int(start0), int(end0), float(value)))

        if len(observed) != len(expected):
            return fail(f"interval count differs: observed={len(observed)} expected={len(expected)}")

        for index, (obs, exp) in enumerate(zip(observed, expected), start=1):
            if obs[:3] != exp[:3] or abs(obs[3] - exp[3]) > tolerance:
                return fail(f"interval {index} differs: observed={obs!r} expected={exp!r}")
    finally:
        bw.close()

    print(f"[pyBigWig compatibility] PASS: {bigwig_path}")
    return 0


def write_file(output_path, bedgraph_path, chrom_sizes_path):
    pyBigWig = import_pybigwig()
    signal = read_bedgraph(bedgraph_path)
    chrom_sizes = read_chrom_sizes(chrom_sizes_path)
    chrom_rank = {chrom: index for index, (chrom, _) in enumerate(chrom_sizes)}
    signal.sort(key=lambda row: (chrom_rank[row[0]], row[1], row[2]))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    bw = pyBigWig.open(str(output_path), "w")
    if bw is None:
        return fail("failed to create output BigWig")
    try:
        bw.addHeader(chrom_sizes)
        if signal:
            bw.addEntries(
                [row[0] for row in signal],
                [row[1] for row in signal],
                ends=[row[2] for row in signal],
                values=[row[3] for row in signal],
            )
    finally:
        bw.close()

    print(f"[pyBigWig compatibility] WROTE: {output_path}")
    return 0


def main(argv):
    if len(argv) < 2:
        return fail("usage: pybigwig-bridge.001.py <check|validate|write> ...")
    command = argv[1]
    try:
        if command == "check":
            pyBigWig = import_pybigwig()
            print(getattr(pyBigWig, "__version__", "unknown"))
            return 0
        if command == "validate" and len(argv) == 5:
            return validate_file(argv[2], argv[3], argv[4])
        if command == "write" and len(argv) == 5:
            return write_file(argv[2], argv[3], argv[4])
        return fail("invalid command or argument count")
    except Exception as exc:
        return fail(str(exc))


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
