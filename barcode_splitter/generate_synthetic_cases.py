# -*- coding: utf-8 -*-
import argparse, gzip, os, random, string
from pathlib import Path

DNA = "ACGT"

def rand_barcode(length):
    return "".join(random.choice(DNA) for _ in range(length))

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

def make_barcode_set(k, length, min_hd=2, max_tries=10000):
    """Generate k barcodes of given length with >= min_hd pairwise Hamming distance."""
    bcs = []
    tries = 0
    while len(bcs) < k and tries < max_tries:
        tries += 1
        cand = rand_barcode(length)
        if any(hamming(cand, x) < min_hd for x in bcs):
            continue
        bcs.append(cand)
    if len(bcs) < k:
        raise RuntimeError("Could not generate barcode set with desired constraints.")
    return bcs

def mutate_barcode(bc, max_mismatches, max_deletions, p_mismatch, p_deletion):
    """Return (mutated_prefix, used_mm, used_del). Deletions remove chars from the barcode (left-aligned)."""
    seq = list(bc)
    used_mm = 0
    used_del = 0

    # Apply deletions first (remove positions)
    for i in range(len(seq)):
        if used_del >= max_deletions:
            break
        if random.random() < p_deletion:
            # delete this char (simulate drop from barcode)
            seq[i] = None
            used_del += 1
    seq = [x for x in seq if x is not None]

    # Apply mismatches
    for i, c in enumerate(seq):
        if used_mm >= max_mismatches:
            break
        if random.random() < p_mismatch:
            choices = [x for x in DNA if x != c]
            seq[i] = random.choice(choices)
            used_mm += 1

    return "".join(seq), used_mm, used_del

def rand_qual(n, qchar="I"):
    # high quality default 'I'; keep simple and fast
    return qchar * n

def write_fastq_gz(path, records):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        for h, s, q in records:
            f.write(h + "\n")
            f.write(s + "\n")
            f.write("+\n")
            f.write(q + "\n")

def generate_case(out_dir: Path, case_id: str, n_reads: int, n_samples: int,
                  bc_len: int, body_len: int,
                  mm_rate: float, del_rate: float, noise_rate: float,
                  max_allowed_mm: int, max_allowed_del: int):
    """Create one case directory with barcodes.txt, reads.fastq.gz, truth.csv."""
    d = out_dir / case_id
    d.mkdir(parents=True, exist_ok=True)

    # 1) barcodes
    bcs = make_barcode_set(n_samples, bc_len, min_hd=max(2, max_allowed_mm+1))
    sample_names = [f"Sample{i+1}" for i in range(n_samples)]
    with open(d / "barcodes.txt", "w", encoding="utf-8") as f:
        for name, bc in zip(sample_names, bcs):
            f.write(f"{name}\t{bc}\n")

    # 2) reads + truth
    records = []
    truth_rows = ["read_id,true_sample"]
    for ridx in range(1, n_reads + 1):
        # select target: either a real sample or "unassigned" (noise)
        if random.random() < noise_rate:
            # no valid barcode at 5' → should be unassigned
            true_label = "unassigned"
            # add random junk prefix length 0..(bc_len-1) that doesn't equal any barcode prefix
            junk_len = random.randint(0, max(0, bc_len - 1))
            junk = "".join(random.choice(DNA) for _ in range(junk_len))
            body = "".join(random.choice(DNA) for _ in range(body_len))
            seq = junk + body
        else:
            idx = random.randrange(n_samples)
            true_label = sample_names[idx]
            bc = bcs[idx]

            # mutate barcode within or beyond allowed limits (mostly within, to be realistic)
            # We bias to produce mostly solvable cases:
            want_mm = min(max_allowed_mm, int(random.random() < 0.5))  # 0 or 1 if allowed>=1
            want_del = min(max_allowed_del, int(random.random() < 0.3))  # 0 or 1 if allowed>=1

            # Then apply per-base randomness around those wants
            mutated, used_mm, used_del = mutate_barcode(
                bc,
                max_mismatches=want_mm,
                max_deletions=want_del,
                p_mismatch=mm_rate,
                p_deletion=del_rate
            )

            body = "".join(random.choice(DNA) for _ in range(body_len))
            seq = mutated + body

        read_id = f"{case_id}_read{ridx:06d}"
        header = f"@{read_id} sample={true_label}"
        qual = rand_qual(len(seq))
        records.append((header, seq, qual))
        truth_rows.append(f"{read_id},{true_label}")

    write_fastq_gz(d / "reads.fastq.gz", records)
    with open(d / "truth.csv", "w", encoding="utf-8") as f:
        f.write("\n".join(truth_rows))

def main():
    ap = argparse.ArgumentParser(description="Generate synthetic demux cases.")
    ap.add_argument("--root", default=".", help="Project root; cases/ will be created here")
    ap.add_argument("--cases", type=int, default=1000, help="Number of cases to generate")
    ap.add_argument("--reads-per-case", type=int, default=300, help="Reads per case")
    ap.add_argument("--samples-per-case", type=int, default=6, help="Barcodes / samples per case")
    ap.add_argument("--barcode-len", type=int, default=6, help="Barcode length")
    ap.add_argument("--read-body-len", type=int, default=80, help="Read body length excluding barcode")
    ap.add_argument("--mismatch-rate", type=float, default=0.05, help="Per-base mismatch rate in barcode")
    ap.add_argument("--deletion-rate", type=float, default=0.00, help="Per-base deletion rate in barcode")
    ap.add_argument("--noise-rate", type=float, default=0.10, help="Fraction of reads with no valid barcode (-> unassigned)")
    ap.add_argument("--max-allowed-mm", type=int, default=1, help="Target: typical allowed mismatches")
    ap.add_argument("--max-allowed-del", type=int, default=0, help="Target: typical allowed deletions")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    args = ap.parse_args()

    random.seed(args.seed)
    root = Path(args.root).resolve()
    cases_dir = root / "cases"
    cases_dir.mkdir(parents=True, exist_ok=True)

    for i in range(1, args.cases + 1):
        case_id = f"case{i:04d}"
        generate_case(
            cases_dir, case_id,
            n_reads=args.reads_per_case,
            n_samples=args.samples_per_case,
            bc_len=args.barcode_len,
            body_len=args.read_body_len,
            mm_rate=args.mismatch_rate,
            del_rate=args.deletion-rate if hasattr(args, "deletion-rate") else args.deletion_rate,  # guard
            noise_rate=args.noise_rate,
            max_allowed_mm=args.max_allowed_mm,
            max_allowed_del=args.max_allowed_del
        )
    print(f"Done. Generated {args.cases} cases under {cases_dir}")

if __name__ == "__main__":
    main()
'''
python ./generate_synthetic_cases.py \
  --root . \
  --cases 1000 \
  --reads-per-case 2000 \
  --samples-per-case 3 \
  --barcode-len 6 \
  --read-body-len 80 \
  --mismatch-rate 0.0 \
  --deletion-rate 0.0 \
  --noise-rate 0.02 \
  --max-allowed-mm 0 \
  --max-allowed-del 0 \
  --seed 1
'''