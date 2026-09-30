#!/usr/bin/env python3
"""
Estimate SLURM resources for `bowtie2-build --large-index` from the size of
the .fna files in a folder.

All sizes are binary units (GiB), matching seff and `du -h`.

Built-in calibration: four completed Dardel jobs (256 cores, 1.72 TB node).

    job                 input GiB   peak mem GiB   wall h   index GiB   contigs
    SE_Plants               67.6          686        10.32      120     2,191,265
    Greenland_Masked       129.0         1444        23.41      249     1,667,698
    NEU_Herps              147.0         1567        28.39      277     1,412,063
    PhyloNorway            151.9         1577        26.97      324   311,635,022

Model
    memory = k * input                (through-origin fit; ~10.7 GiB per GiB)
    time   = c * input ** p           (power-law fit, p ~ 1.24, when >= 3 jobs)
             or a linear fit through the origin with fewer points
    index  = a * input + b * contigs_millions   (no-intercept 2-variable fit;
             a ~ 1.89 GiB/GiB, b ~ 0.12 GiB/million contigs. PhyloNorway's huge
             contig count, from many short sequences, is what pins down b;
             falls back to a flat ~1.98 GiB/GiB ratio when contig counts aren't
             available for calibration or for the target folder)

Add your own jobs with --calibration jobs.json:
    [{"name": "x", "input_gb": 147, "mem_gb": 1567, "wall_hours": 28.39,
      "cpu_hours": 725.2, "cores": 256, "index_gb": 277, "contigs": 1412063}, ...]

Two ways to run this:
  1. CLI:
       python est_index_resources.py /path/to/fna_dir --emit-sbatch
  2. Snakemake `script:` directive (see bottom of file) — reads every path in
     snakemake.input as an individual .fna file (e.g. from an `expand(...)`
     list) and writes a TSV report to snakemake.output[0].
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

FNA_SUFFIXES = (".fna", ".fa", ".fasta")
GIB = 1024 ** 3

# Defaults matching the Dardel "memory" partition profile: whole node granted
# regardless of requested cpus, ~1760 GiB memory, 5-day (7200 min) time limit.
DEFAULT_NODE_MEM_GIB = 1760.0
DEFAULT_MAX_HOURS = 120.0
DEFAULT_PARTITION = "memory"


def h(d, hh, m, s):
    return d * 24 + hh + m / 60 + s / 3600


DEFAULT_JOBS = [
    {"name": "SE_Plants", "input_gb": 67.6, "mem_gb": 686.42, "index_gb": 120.0, "contigs": 2191265,
     "wall_hours": h(0, 10, 19, 26), "cpu_hours": h(12, 16, 56, 20), "cores": 256},
    {"name": "Greenland_Masked", "input_gb": 129.0, "mem_gb": 1.41 * 1024, "index_gb": 249.0, "contigs": 1667698,
     "wall_hours": h(0, 23, 24, 18), "cpu_hours": h(27, 0, 16, 28), "cores": 256},
    {"name": "NEU_Herps", "input_gb": 147.0, "mem_gb": 1.53 * 1024, "index_gb": 277.0, "contigs": 1412063,
     "wall_hours": h(1, 4, 23, 29), "cpu_hours": h(30, 5, 10, 23), "cores": 256},
    {"name": "PhyloNorway", "input_gb": 151.9, "mem_gb": 1.54 * 1024, "index_gb": 324.0, "contigs": 311635022,
     "wall_hours": h(1, 2, 58, 14), "cpu_hours": h(29, 8, 32, 41), "cores": 256},
]


# --------------------------------------------------------------------------
# Scanning helpers
# --------------------------------------------------------------------------
# `scan`, `count_sequences`, and `count_contigs` all take an explicit list of
# .fna/.fa/.fasta file paths — the caller resolves what those paths are,
# whether that's every file matching a glob under a folder (CLI/`--recursive`)
# or an already-enumerated list Snakemake handed in via `input:`.

def resolve_folder_files(folder: Path, recursive: bool) -> list[Path]:
    """CLI helper: find every .fna/.fa/.fasta (optionally .gz) file under a folder."""
    it = folder.rglob("*") if recursive else folder.glob("*")
    files = []
    for p in it:
        if not p.is_file():
            continue
        name = p.name.lower()
        base = name[:-3] if name.endswith(".gz") else name
        if base.endswith(FNA_SUFFIXES):
            files.append(p)
    return files


def scan(files: list[Path], gz_factor: float):
    """Return (n_files, on_disk_bytes, estimated_uncompressed_bytes)."""
    n, disk, est = 0, 0, 0.0
    for p in files:
        size = p.stat().st_size
        gz = p.name.lower().endswith(".gz")
        n += 1
        disk += size
        est += size * gz_factor if gz else size
    return n, disk, est


def count_sequences(files: list[Path]):
    """Optional slow pass: count '>' records and bases (informational only)."""
    import gzip
    seqs = bases = 0
    for p in files:
        opener = gzip.open if p.name.lower().endswith(".gz") else open
        with opener(p, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    seqs += 1
                else:
                    bases += len(line.strip())
    return seqs, bases


def count_contigs(files: list[Path]):
    """Fast pass: count '>' header lines only (no base counting)."""
    import gzip
    seqs = 0
    for p in files:
        opener = gzip.open if p.name.lower().endswith(".gz") else open
        with opener(p, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    seqs += 1
    return seqs


# --------------------------------------------------------------------------
# Fitting helpers
# --------------------------------------------------------------------------

def fit_bivariate_no_intercept(x1, x2, y):
    """Least squares y = a*x1 + b*x2 (no intercept). Returns (a, b)."""
    s11 = sum(v * v for v in x1)
    s22 = sum(v * v for v in x2)
    s12 = sum(a * b for a, b in zip(x1, x2))
    s1y = sum(a * b for a, b in zip(x1, y))
    s2y = sum(a * b for a, b in zip(x2, y))
    det = s11 * s22 - s12 * s12
    if abs(det) < 1e-9:
        return fit_through_origin(x1, y), 0.0
    a = (s1y * s22 - s2y * s12) / det
    b = (s11 * s2y - s12 * s1y) / det
    return a, b


def fit_through_origin(xs, ys):
    return sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)


def fit_linear(xs, ys):
    """Least-squares y = a*x + b. Returns (a, b)."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    a = num / den if den else 0.0
    return a, my - a * mx


def fit_power_law(xs, ys):
    """Least-squares fit of y = c * x**p in log space. Returns (c, p)."""
    lx, ly = [math.log(x) for x in xs], [math.log(y) for y in ys]
    n = len(xs)
    mx, my = sum(lx) / n, sum(ly) / n
    p = sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx)
    return math.exp(my - p * mx), p


def fmt_time(hours: float) -> str:
    total = int(math.ceil(hours * 3600))
    d, rem = divmod(total, 86400)
    hh, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    return f"{d}-{hh:02d}:{m:02d}:{s:02d}"


# --------------------------------------------------------------------------
# Core estimation (shared by the CLI and the Snakemake entry point)
# --------------------------------------------------------------------------

def estimate(
    files: list[Path],
    *,
    label: str | None = None,
    gz_factor: float = 4.0,
    calibration_path: Path | None = None,
    mem_safety: float = 1.10,
    time_safety: float = 3.0,
    cpus: int | None = None,
    node_mem_gb: float | None = None,
    max_hours: float | None = None,
    count_seqs: bool = False,
    contigs: int | None = None,
    no_contig_count: bool = False,
    partition: str | None = None,
    whole_node: bool = True,
) -> dict:
    """Run the model fit over an explicit list of .fna/.fa/.fasta files and
    return a plain dict of results. `label` is just a display string for the
    report (e.g. the source folder, or a project name) — pass whatever's
    meaningful for the caller. node_mem_gb/max_hours/partition fall back to
    the Dardel "memory" partition's known values when not given. whole_node=True
    (the default) adds a note that --cpus doesn't need tuning on a partition
    that grants a whole node regardless of what's requested.
    Raises SystemExit if `files` is empty."""
    if not files:
        raise SystemExit("No .fna/.fa/.fasta files given.")
    files = [Path(f) for f in files]
    if node_mem_gb is None:
        node_mem_gb = DEFAULT_NODE_MEM_GIB
    if max_hours is None:
        max_hours = DEFAULT_MAX_HOURS
    if partition is None:
        partition = DEFAULT_PARTITION

    jobs = json.loads(calibration_path.read_text()) if calibration_path else DEFAULT_JOBS
    xs = [j["input_gb"] for j in jobs]
    mem_per_gb = fit_through_origin(xs, [j["mem_gb"] for j in jobs])

    idx_jobs = [j for j in jobs if "index_gb" in j]
    idx_per_gb = fit_through_origin([j["input_gb"] for j in idx_jobs], [j["index_gb"] for j in idx_jobs]) \
        if idx_jobs else None
    idx_contig_jobs = [j for j in idx_jobs if "contigs" in j]
    idx_a = idx_d = None
    if len(idx_contig_jobs) >= 2:
        idx_a, idx_d = fit_bivariate_no_intercept(
            [j["input_gb"] for j in idx_contig_jobs],
            [j["contigs"] / 1e6 for j in idx_contig_jobs],
            [j["index_gb"] for j in idx_contig_jobs],
        )

    if len(set(xs)) >= 3:
        t_c, t_p = fit_power_law(xs, [j["wall_hours"] for j in jobs])
        time_desc = f"{t_c:.4f} * GiB^{t_p:.2f} h"
    else:
        t_c, t_p = fit_through_origin(xs, [j["wall_hours"] for j in jobs]), 1.0
        time_desc = f"{t_c:.3f} h/GiB (linear)"

    avg_cores = [j["cpu_hours"] / j["wall_hours"] for j in jobs if "cpu_hours" in j]
    if len(avg_cores) == len(jobs) and len(jobs) >= 2:
        core_a, core_b = fit_linear(xs, avg_cores)
        core_desc = f"{core_a:+.3f}*GiB {'+' if core_b >= 0 else '-'} {abs(core_b):.1f} cores"
    else:
        core_a, core_b = 0.0, (sum(avg_cores) / len(avg_cores) if avg_cores else 32.0)
        core_desc = f"{core_b:.1f} cores (flat, insufficient data to fit a trend)"

    n, disk, est = scan(files, gz_factor)
    if n == 0:
        raise SystemExit("No .fna/.fa/.fasta files found.")
    input_gb = est / GIB

    mem_point = input_gb * mem_per_gb
    mem_gb = mem_point * mem_safety
    hours = t_c * input_gb ** t_p * time_safety

    est_cpus = core_a * input_gb + core_b
    req_cpus = cpus if cpus is not None else max(1, min(256, round(est_cpus)))

    n_seqs = n_bases = None
    if count_seqs:
        n_seqs, n_bases = count_sequences(files)

    n_contigs = contigs
    index_gb = None
    index_model_desc = None
    if idx_a is not None:
        if n_contigs is None and not no_contig_count:
            n_contigs = count_contigs(files)
        if n_contigs is not None:
            index_gb = idx_a * input_gb + idx_d * (n_contigs / 1e6)
            index_model_desc = f"{idx_a:.2f} GiB/GiB input + {idx_d:.3f} GiB/M contigs"
        elif idx_per_gb:
            index_gb = input_gb * idx_per_gb
            index_model_desc = f"{idx_per_gb:.2f} GiB/GiB input (contig count unavailable)"
    elif idx_per_gb:
        index_gb = input_gb * idx_per_gb
        index_model_desc = f"{idx_per_gb:.2f} GiB/GiB input"

    # Notes are scoped to the report section they belong under, so the CLI
    # can print each one right after the line it explains.
    notes = []  # list of (section, level, text)
    lo, hi = min(xs), max(xs)
    if input_gb < 0.5 * lo or input_gb > 1.1 * hi:
        notes.append(("input", "NOTE", f"{input_gb:.0f} GiB is outside the calibrated range "
                                        f"({lo:.0f}-{hi:.0f} GiB); the estimate is an extrapolation."))
    mem_capped = False
    if mem_point > node_mem_gb:
        max_in = node_mem_gb / mem_per_gb
        notes.append(("mem", "WARNING", f"needs ~{mem_point:.0f} GiB before any margin but a node has "
                                         f"{node_mem_gb:.0f} GiB (largest input at this fit: ~{max_in:.0f} GiB). "
                                         f"Split the reference set, or trade speed for memory "
                                         f"(bowtie2-build --bmax / --dcv)."))
        mem_capped = True
    elif mem_gb > node_mem_gb:
        notes.append(("mem", "NOTE", f"the estimate with margin ({mem_gb:.0f} GiB) exceeds the node "
                                      f"({node_mem_gb:.0f} GiB); --mem is capped at the node size, "
                                      f"headroom is thin."))
        mem_capped = True
    if max_hours and hours > max_hours:
        notes.append(("time", "WARNING", f"estimated {hours:.1f} h exceeds the partition limit of "
                                          f"{max_hours} h."))
    if req_cpus < est_cpus - 1:
        notes.append(("cpus", "NOTE", f"requesting {req_cpus} cores but bowtie2-build looks to actually use "
                                       f"about {est_cpus:.0f} at this input size; fewer cores may slow the "
                                       f"job beyond the time estimate, which assumes cores aren't the "
                                       f"bottleneck."))
    if whole_node:
        notes.append(("cpus", "NOTE", "there is no need to set the number of cpus if your HPC grants an "
                                       "entire node on a memory partition."))
    if index_gb is not None:
        notes.append(("index", "INFO", "Check free disk space at the output path before running "
                                        "bowtie2-build --large-index."))

    if label is None:
        try:
            import os
            label = os.path.commonpath([str(f) for f in files]) if len(files) > 1 else str(files[0])
        except ValueError:
            label = f"{len(files)} input files"

    return {
        "folder": label,
        "n_files": n,
        "on_disk_gib": disk / GIB,
        "input_gib": input_gb,
        "n_sequences": n_seqs,
        "n_bases": n_bases,
        "n_contigs": n_contigs,
        "calibration_n_jobs": len(jobs),
        "mem_per_gib": mem_per_gb,
        "time_model": time_desc,
        "core_model": core_desc,
        "idx_per_gib": idx_per_gb,
        "index_model_desc": index_model_desc,
        "mem_point_gib": mem_point,
        "mem_gib": min(mem_gb, node_mem_gb),
        "mem_mb_slurm": math.ceil(min(mem_gb, node_mem_gb) * 1074),
        "mem_capped": mem_capped,
        "time_hours": hours,
        "time_slurm": fmt_time(hours),
        "cpus": req_cpus,
        "est_cpus": est_cpus,
        "index_gib": index_gb,
        "partition": partition,
        "notes": notes,
    }


# --------------------------------------------------------------------------
# Human-readable report + sbatch header (CLI use)
# --------------------------------------------------------------------------

def _print_notes(r, section):
    for sec, level, text in r["notes"]:
        if sec == section:
            print(f"{level}: {text}")


def print_report(r: dict, emit_sbatch: bool = False, partition: str | None = None):
    partition = partition or r.get("partition")

    print(f"Files found          : {r['n_files']}")
    print(f"On-disk size         : {r['on_disk_gib']:.1f} GiB")
    print(f"Est. uncompressed    : {r['input_gib']:.1f} GiB")
    _print_notes(r, "input")
    if r["n_sequences"] is not None:
        print(f"Sequences / bases    : {r['n_sequences']:,} / {r['n_bases']:,}")

    idx_part = f"; index {r['idx_per_gib']:.2f} GiB/GiB input" if r["idx_per_gib"] else ""
    print(f"Calibration ({r['calibration_n_jobs']} jobs): mem {r['mem_per_gib']:.2f} GiB/GiB input; "
          f"time {r['time_model']}; cores {r['core_model']}{idx_part}")

    print(f"Suggested --mem      : {r['mem_mb_slurm']}M  (~{r['mem_gib']:.0f} GiB)")
    _print_notes(r, "mem")
    print(f"Suggested --time     : {r['time_slurm']}")
    _print_notes(r, "time")
    print(f"Suggested --cpus     : {r['cpus']}  (fitted avg-cores-used estimate: {r['est_cpus']:.1f})")
    _print_notes(r, "cpus")

    if r["index_gib"] is not None:
        if r["n_contigs"] is not None:
            print(f"Contigs (headers)    : {r['n_contigs']:,}")
        index_notes = "; ".join(text for sec, _, text in r["notes"] if sec == "index")
        suffix = f"; {index_notes}" if index_notes else ""
        print(f"Est. index size      : ~{r['index_gib']:.0f} GiB ({r['index_model_desc']}{suffix})")

    if emit_sbatch:
        lines = [
            "#!/bin/bash",
            "#SBATCH --job-name=bt2build",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            f"#SBATCH --cpus-per-task={r['cpus']}",
            f"#SBATCH --mem={r['mem_mb_slurm']}M",
            f"#SBATCH --time={r['time_slurm']}",
        ]
        if partition:
            lines.append(f"#SBATCH --partition={partition}")
        print()
        print("\n".join(lines))


def write_tsv(r: dict, output_path: Path):
    """Write the estimate as a two-column (metric, value) TSV report."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    all_notes = "; ".join(f"[{sec}/{level}] {text}" for sec, level, text in r["notes"])
    rows = [
        ("Folder=", r["folder"]),
        ("Number of files=", r["n_files"]),
        ("Input .fna size (GiB)=", f"{r['on_disk_gib']:.2f}"),
        ("Uncompressed .fna size (GiB)=", f"{r['input_gib']:.2f}"),
        ("Number of contigs=", r["n_contigs"] if r["n_contigs"] is not None else ""),
        (),
        ("Time Est. model=", r["time_model"]),
        ("Core Est. model=", r["core_model"]),
        ("Index Est. model=", r["index_model_desc"] or ""),
        (),
        ("SNAKEMAKE SETTINGS (profiles/config.yaml)",),
        ("Suggested mem_mb=", f"{r['mem_gib']*1074:.1f}"),
        ("Suggested runtime (minutes)=", f"{r['time_hours'] * 60:.2f}"),
        ("Suggested cpus_per_task=", r["cpus"]),
        (),
        ("SLURM SETTINGS EQUIVALENT",),
        ("Suggested --mem (Mb)=", f"{r['mem_mb_slurm']}"),
        ("Suggested --time=", r["time_slurm"]),
        ("Suggested --cpus-per-task=", r["cpus"]),
        ("NOTE: there is no need to set the number of cpus if your HPC grants an entire node on a memory partition.",),
        (),
        ("Est. Index Size (GiB)=", f"{r['index_gib']:.1f}" if r["index_gib"] is not None else ""),
        ("NOTE: Check free disk space at the output path before running bowtie2 build.",),
        (),
        ("", all_notes),
    ]
    with open(output_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(rows)


# --------------------------------------------------------------------------
# CLI entry point
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", type=Path)
    ap.add_argument("--recursive", action="store_true")
    ap.add_argument("--gz-factor", type=float, default=4.0,
                    help="assumed uncompressed/compressed ratio for .gz inputs (default 4.0)")
    ap.add_argument("--calibration", type=Path, help="JSON list of observed jobs (replaces the built-in four)")
    ap.add_argument("--mem-safety", type=float, default=1.10,
                    help="memory multiplier over the fitted peak (default 1.10; worst observed error is ~5%%)")
    ap.add_argument("--time-safety", type=float, default=1.25,
                    help="wall-time multiplier over the fitted time (default 1.25)")
    ap.add_argument("--cpus", type=int, default=None,
                    help="cores to request (default: estimated from input size, since bowtie2-build "
                         "only kept 10-12%% of 256 cores busy on the calibration jobs)")
    ap.add_argument("--node-mem-gb", type=float, default=None,
                    help=f"max memory per node in GiB (default {DEFAULT_NODE_MEM_GIB:.0f})")
    ap.add_argument("--max-hours", type=float, default=None,
                    help=f"partition wall-time limit in hours (default {DEFAULT_MAX_HOURS:.0f})")
    ap.add_argument("--count-seqs", action="store_true",
                    help="also count sequences/bases (slow; informational only)")
    ap.add_argument("--contigs", type=int, default=None,
                    help="known contig/sequence count for the folder, to feed the index-size model "
                         "(default: counted automatically unless --no-contig-count is given)")
    ap.add_argument("--no-contig-count", action="store_true",
                    help="skip scanning for contig count; index size then falls back to the "
                         "input-size-only model")
    ap.add_argument("--no-whole-node-note", dest="whole_node", action="store_false",
                    help="suppress the note saying --cpus doesn't need tuning on a whole-node partition")
    ap.add_argument("--emit-sbatch", action="store_true")
    ap.add_argument("--partition", default=None, help=f"default: {DEFAULT_PARTITION!r}")
    ap.add_argument("--output", type=Path, default=None,
                    help="also write the estimate as a TSV report to this path")
    args = ap.parse_args()

    files = resolve_folder_files(args.folder, args.recursive)
    if not files:
        raise SystemExit(f"No .fna/.fa/.fasta files found under {args.folder}")

    r = estimate(
        files,
        label=str(args.folder),
        gz_factor=args.gz_factor,
        calibration_path=args.calibration,
        mem_safety=args.mem_safety,
        time_safety=args.time_safety,
        cpus=args.cpus,
        node_mem_gb=args.node_mem_gb,
        max_hours=args.max_hours,
        count_seqs=args.count_seqs,
        contigs=args.contigs,
        no_contig_count=args.no_contig_count,
        partition=args.partition,
        whole_node=args.whole_node,
    )
    print_report(r, emit_sbatch=args.emit_sbatch, partition=r["partition"])
    if args.output:
        write_tsv(r, args.output)

def run_from_snakemake(snakemake) -> None:
    input_files = [Path(f).expanduser() for f in snakemake.input]
    output_path = Path(snakemake.output[0]).expanduser()

    missing = [f for f in input_files if not f.exists()]
    if missing:
        raise SystemExit(f"Input file(s) not found: {', '.join(str(f) for f in missing)}")

    params = getattr(snakemake, "params", {})

    def p(name, default):
        try:
            return params[name]
        except (KeyError, TypeError, AttributeError):
            return getattr(params, name, default)

    calibration = p("calibration", None)
    r = estimate(
        input_files,
        label=p("label", None),
        gz_factor=p("gz_factor", 4.0),
        calibration_path=Path(calibration) if calibration else None,
        mem_safety=p("mem_safety", 1.10),
        time_safety=p("time_safety", 3.0),
        cpus=p("cpus", None),
        node_mem_gb=p("node_mem_gb", None),
        max_hours=p("max_hours", None),
        count_seqs=False,
        contigs=p("contigs", None),
        no_contig_count=p("no_contig_count", False),
        partition=p("partition", None),
        whole_node=p("whole_node", True),
    )
    write_tsv(r, output_path)


try:
    snakemake  # noqa: F821 — injected by Snakemake's `script:` directive
except NameError:
    snakemake = None

if snakemake is not None:
    run_from_snakemake(snakemake)
elif __name__ == "__main__":
    main()