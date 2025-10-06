

import argparse
import gzip
import os
from collections import namedtuple, defaultdict
from typing import List, Tuple, Optional, Dict, Iterable

Barcode = namedtuple("Barcode", ["name", "seq", "len"])

def read_barcode_table(path: str) -> List[Barcode]:
    out = []
    with open(path, "r", encoding="utf-8") as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            name, seq = parts[0], parts[1].upper().replace("U", "T")
            out.append(Barcode(name=name, seq=seq, len=len(seq)))
    if not out:
        raise ValueError("No barcodes parsed from file.")
    return out

def open_textstream(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")

def write_fastq_record(handle, h, s, p, q):
    handle.write(h); handle.write(s); handle.write(p); handle.write(q)

def parse_fastq(stream) -> Iterable[Tuple[str, str, str, str]]:
    while True:
        h = stream.readline()
        if not h:
            break
        s = stream.readline()
        p = stream.readline()
        q = stream.readline()
        if not (s and p and q):
            raise ValueError("Truncated FASTQ input.")
        yield h, s, p, q

def match_barcode_at_5prime(read_seq: str, barcode: str, max_mismatches: int, max_deletions: int) -> Optional[Tuple[int, int]]:
    B = len(barcode)
    R = len(read_seq)
    INF = (10**9, 10**9)
    DP = [[INF]*(R+1) for _ in range(B+1)]
    DP[0][0] = (0, 0)
    for i in range(B):
        for j in range(R+1):
            mm, dd = DP[i][j]
            if mm == 10**9:
                continue
            if dd + 1 <= max_deletions:
                if DP[i+1][j] > (mm, dd+1):
                    DP[i+1][j] = (mm, dd+1)
            if j < R:
                add_mis = 0 if barcode[i] == read_seq[j] else 1
                if mm + add_mis <= max_mismatches:
                    cand = (mm + add_mis, dd)
                    if DP[i+1][j+1] > cand:
                        DP[i+1][j+1] = cand
    best = None
    for j in range(R+1):
        mm, dd = DP[B][j]
        if mm <= max_mismatches and dd <= max_deletions:
            if best is None or (mm, dd, j) < (best[1], best[2], best[0]):
                best = (j, mm, dd)
    if best is None:
        return None
    used_read_len, mm, _dd = best
    return used_read_len, mm

def find_best_barcode(read_seq: str, barcodes: List[Barcode], max_mismatches: int, max_deletions: int) -> Optional[Tuple[Barcode, int]]:
    best = None
    for idx, bc in enumerate(barcodes):
        read_prefix = read_seq[:bc.len]
        hit = match_barcode_at_5prime(read_prefix, bc.seq, max_mismatches, max_deletions)
        if hit is None:
            continue
        consumed_len, mm = hit
        dd = bc.len - consumed_len
        score = (mm, dd, -consumed_len, idx)
        if best is None or score < best[0]:
            best = (score, bc, consumed_len)
    if best is None:
        return None
    _, bc, consumed_len = best
    return bc, consumed_len

def run(fastq_path: str, barcode_path: str, outdir: str, mismatches: int, deletions: int):
    os.makedirs(outdir, exist_ok=True)
    barcodes = read_barcode_table(barcode_path)
    writers: Dict[str, any] = {}
    def _open_out(name: str):
        fn = os.path.join(outdir, f"{name}.fastq.gz")
        fh = gzip.open(fn, "wt")
        writers[name] = fh
        return fh
    for bc in barcodes:
        _open_out(bc.name)
    unassigned = _open_out("unassigned")
    counts = defaultdict(int)
    with open_textstream(fastq_path, "rt") as fh:
        for h, s, p, q in parse_fastq(fh):
            seq = s.strip().upper().replace("U", "T")
            qual = q.strip()
            hit = find_best_barcode(seq, barcodes, mismatches, deletions)
            if hit is None:
                write_fastq_record(unassigned, h, s, p, q)
                counts["unassigned"] += 1
                continue
            bc, consumed = hit
            trimmed_seq = seq[consumed:]
            trimmed_qual = qual[consumed:]
            write_fastq_record(writers[bc.name], h, trimmed_seq + "\n", p, trimmed_qual + "\n")
            counts[bc.name] += 1
    for fh in writers.values():
        fh.close()
    total = sum(v for k, v in counts.items())
    lines = ["# Barcode Splitter summary"]
    for bc in barcodes:
        lines.append(f"{bc.name}\t{counts[bc.name]}")
    lines.append(f"unassigned\t{counts['unassigned']}")
    lines.append(f"total\t{total}")
    with open(os.path.join(outdir, "summary.tsv"), "w", encoding="utf-8") as sf:
        sf.write("\n".join(lines) + "\n")

def main():
    ap = argparse.ArgumentParser(description="Split a FASTQ by 5' barcodes with allowed mismatches and barcode deletions.")
    ap.add_argument("--fastq", required=True, help="Input FASTQ(.gz)")
    ap.add_argument("--barcodes", required=True, help="Barcode table (id\\tseq)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--mismatches", type=int, default=0, help="Allowed mismatches")
    ap.add_argument("--deletions", type=int, default=0, help="Allowed barcode nucleotide deletions")
    args = ap.parse_args()
    run(fastq_path=args.fastq, barcode_path=args.barcodes, outdir=args.outdir, mismatches=args.mismatches, deletions=args.deletions)

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