import argparse, csv, gzip, os, re, sys, json, shutil, itertools
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
from typing import Optional
# -*- coding: utf-8 -*-

def discover_cases(root: Path, cases_dir: str = "cases"):
    """Discover cases under root/cases: expect case*/reads.fastq(.gz) + barcodes.txt (+ truth.csv optional)"""
    rows = []
    base = root / cases_dir
    if not base.exists():
        return rows
    for d in sorted(base.glob("case*")):
        if not d.is_dir():
            continue
        fq = None
        for cand in ["reads.fastq.gz", "reads.fastq"]:
            p = d / cand
            if p.exists():
                fq = p
                break
        bc = d / "barcodes.txt"
        tr = d / "truth.csv"
        if fq and bc.exists():
            rows.append({
                "case_id": d.name,
                "fastq": str(fq.relative_to(root)),
                "barcodes": str(bc.relative_to(root)),
                "truth_csv": str(tr.relative_to(root)) if tr.exists() else ""
            })
    return rows

def write_manifest(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["case_id","fastq","barcodes","truth_csv"])
        w.writeheader(); w.writerows(rows)


def read_summary_counts(outdir: Path):
    """Return (counts_dict, total) by preferring summary.tsv; fallback to counting FASTQ records."""
    import csv, glob, gzip
    sp = outdir / "summary.tsv"
    counts = {}
    total = None
    # 1) Try summary.tsv
    if sp.exists():
        with open(sp, "r", encoding="utf-8", errors="ignore") as f:
            head = f.read(4096)
            delim = "	" if ("	" in head) else ("," if ("," in head) else "	")
            f.seek(0)
            rd = csv.DictReader(f, delimiter=delim)
            for row in rd:
                key = (row.get("sample") or row.get("name") or row.get("barcode") or "").strip() or "unknown"
                try:
                    c = int(row.get("count", 0))
                except Exception:
                    c = 0
                counts[key] = counts.get(key, 0) + c
        total = sum(counts.values())
        if total and total > 0:
            return counts, total
    # 2) Fallback: count reads from *.fastq / *.fastq.gz (lines/4)
    files = glob.glob(str(outdir / "*.fastq")) + glob.glob(str(outdir / "*.fastq.gz"))
    for fq in files:
        name = Path(fq).name
        sample = name.split(".fastq")[0]
        sample = "unassigned" if sample.lower() == "unassigned" else sample
        try:
            if fq.endswith(".gz"):
                with gzip.open(fq, "rt", encoding="utf-8", errors="ignore") as fh:
                    nlines = sum(1 for _ in fh)
            else:
                with open(fq, "rt", encoding="utf-8", errors="ignore") as fh:
                    nlines = sum(1 for _ in fh)
            nreads = nlines // 4
            if nreads > 0:
                counts[sample] = counts.get(sample, 0) + nreads
        except Exception:
            pass
    total = sum(counts.values()) if counts else None
    return counts, total

def run_one(case_row, args):
    root = Path(args.root).resolve()
    case_id = case_row["case_id"]
    fq = (root / case_row["fastq"]).resolve()
    bc = (root / case_row["barcodes"]).resolve()
    truth_csv = (root / case_row["truth_csv"]).resolve() if case_row.get("truth_csv") else None

    outdir = Path(args.outdir) / case_id
    outdir.mkdir(parents=True, exist_ok=True)

    # run splitter
    cmd = [
        sys.executable, args.splitter,
        "--fastq", str(fq),
        "--barcodes", str(bc),
        "--outdir", str(outdir),
        "--mismatches", str(args.mismatches),
        "--deletions", str(args.deletions),
    ]
    if args.extra and args.extra.strip():
        cmd.extend(args.extra.strip().split())

    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
        rc = r.returncode
        stderr = (r.stderr or "").strip()
        stdout = (r.stdout or "").strip()
    except Exception as e:
        return {
            "case_id": case_id, "status": "ERROR",
            "retcode": -1, "stderr": str(e)[:600],
            "total_reads": "", "unassigned_reads": "", "summary_json": "{}"
        }

    # parse summary to decide success
    counts, total = read_summary_counts(outdir)
    unassigned_reads = counts.get("unassigned", 0)
    summary_json = json.dumps(counts, ensure_ascii=False)

    # define success: rc==0 AND summary.tsv exists AND total>0
    status = "OK" if (rc == 0 and total and total > 0) else "ERROR"
    # shorten stderr for csv
    if len(stderr) > 600: stderr = stderr[:600]

    return {
        "case_id": case_id,
        "status": status,
        "retcode": rc,
        "stderr": stderr,
        "total_reads": total if total is not None else "",
        "unassigned_reads": unassigned_reads if total is not None else "",
        "summary_json": summary_json
    }

def read_manifest(path: Path):
    rows = []
    if not path.exists():
        return rows
    with open(path, "r", newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for row in rd:
            if row.get("case_id"):
                rows.append(row)
    return rows

def main():
    ap = argparse.ArgumentParser(description="Auto-generate manifest, run cases, and log results to CSV.")
    ap.add_argument("--splitter", default="Barcode Splitter.py", help="Path to splitter script")
    ap.add_argument("--root", default=".", help="Project root; manifest paths are relative to this")
    ap.add_argument("--cases-dir", default="cases", help="Directory under root that contains case* folders")
    ap.add_argument("--manifest", default="manifest.csv", help="Write or reuse this manifest file")
    ap.add_argument("--outdir", default="demux_results", help="Output dir for per-case results")
    ap.add_argument("--results", default="run_results.csv", help="CSV file to write per-case run status")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--mismatches", type=int, default=0)
    ap.add_argument("--deletions", type=int, default=0)
    ap.add_argument("--extra", default="", help="Extra args to pass to splitter")
    ap.add_argument("--reuse-manifest", action="store_true", help="Reuse existing manifest.csv (do not rescan)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    manifest_path = Path(args.manifest).resolve()

    # 1) make or reuse manifest
    if args.reuse_manifest and manifest_path.exists():
        rows = read_manifest(manifest_path)
        print(f"[info] Reusing manifest: {manifest_path} ({len(rows)} cases)")
    else:
        rows = discover_cases(root, args.cases_dir)
        write_manifest(rows, manifest_path)
        print(f"[info] Wrote manifest: {manifest_path} ({len(rows)} cases)")

    if not rows:
        print("[warn] No cases discovered. Check --root and --cases-dir.")
        return

    # 2) run in parallel
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_one, row, args): row["case_id"] for row in rows}
        for fut in as_completed(futs):
            res = fut.result()
            results.append(res)
            cid = res["case_id"]
            print(f"[{cid}] {res['status']} ret={res['retcode']} total={res.get('total_reads','')} unassigned={res.get('unassigned_reads','')}")

    # 3) write CSV
    fieldnames = ["case_id","status","retcode","total_reads","unassigned_reads","summary_json","stderr"]
    with open(args.results, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in sorted(results, key=lambda x: x["case_id"]):
            w.writerow({k: r.get(k, "") for k in fieldnames})

    # 4) print aggregate
    ok = sum(1 for r in results if r["status"] == "OK")
    err = len(results) - ok
    print(f"\n[done] Cases total={len(results)}, success={ok}, failed={err}.")
    print(f"[done] Per-case results -> {Path(args.results).resolve()}")
    print(f"[done] Manifest -> {manifest_path}")
    print(f"[done] Outputs per case under -> {Path(args.outdir).resolve()}")

if __name__ == "__main__":
    main()
"""
python ./auto_run_cases.py \
  --splitter "./Barcode Splitter.py" \
  --root /Workspace10/zehuasun/Biobuild/barcode_splitter \
  --cases-dir cases \
  --outdir demux_results \
  --results run_results.csv \
  --mismatches 1 \
  --deletions 0 \
  --workers 8

"""