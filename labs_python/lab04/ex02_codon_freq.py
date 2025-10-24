#!/usr/bin/env python3
import argparse
from collections import Counter
from pathlib import Path
import json
import gzip
import matplotlib.pyplot as plt

# --- Genetic code (RNA codons -> 1-letter amino acids). Stops as '*'.
GENETIC_CODE = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}
START = "AUG"
STOPS = {"UAA","UAG","UGA"}

# Pretty names for AAs we care about in this lab
AA_NAMES = {
    "L": "Leucine",
    "V": "Valine",
    "T": "Threonine",
    "E": "Glutamic acid (glutamate)",
    "S": "Serine",
}

# Minimal nutrition KB to embed in the output (examples, not exhaustive)
NUTRITION_KB = {
    "L": {
        "low_foods": ["Fruits (e.g., apples)", "Oils/fats (e.g., olive oil)"],
        "rich_foods": ["Chicken breast", "Firm tofu/soy", "Beef", "Fish", "Dairy"],
    },
    "V": {
        "low_foods": ["Fruits", "Oils/fats"],
        "rich_foods": ["Chicken breast", "Beef", "Eggs", "Dairy", "Soy products"],
    },
    "T": {
        "low_foods": ["Fruits", "Oils/fats"],
        "rich_foods": ["Chicken", "Dairy", "Eggs", "Soy"],
    },
    "E": {
        "low_foods": ["Fruits", "Oils/fats"],
        "rich_foods": ["Parmesan/hard cheeses", "Chicken breast", "Soy products"],
    },
    "S": {
        "low_foods": ["Fruits", "Oils/fats"],
        "rich_foods": ["Eggs", "Chicken breast", "Dairy", "Legumes"],
    },
}
NUTRITION_SOURCES = [
    "USDA FoodData Central",
    "MyFoodData (aggregates USDA FDC)",
]

def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def read_fasta_records(path):
    header, seq = None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)

def cat_rna(records):
    s = []
    for _, seq in records:
        s.append(seq.upper().replace("T","U"))
    return "".join(s)

def find_orfs_rna(rna: str, min_len_nt: int = 90):
    orfs = []
    n = len(rna)
    for frame in range(3):
        i = frame
        while i <= n - 3:
            if rna[i:i+3] == START:
                j = i + 3
                while j <= n - 3:
                    c = rna[j:j+3]
                    if c in STOPS:
                        orf = rna[i:j]  # exclude stop
                        if len(orf) >= min_len_nt:
                            orfs.append(orf)
                        i = j + 3
                        break
                    j += 3
                else:
                    # no stop; take to end (trim to multiple of 3)
                    if n - i >= min_len_nt:
                        orfs.append(rna[i: n - (n - i) % 3])
                    i = n
            else:
                i += 3
    return orfs

def count_codons(orfs):
    cnt = Counter()
    for orf in orfs:
        for k in range(0, len(orf) - 2, 3):
            cod = orf[k:k+3]
            if "N" in cod:
                continue
            cnt[cod] += 1
    return cnt

def top_k(cnt, k=10):
    return cnt.most_common(k)

def aa_counts(cnt):
    aac = Counter()
    for cod, v in cnt.items():
        aa = GENETIC_CODE.get(cod, "X")
        if aa not in {"*", "X"}:
            aac[aa] += v
    return aac

def plot_top10(cnt, title, out_path):
    t = top_k(cnt, 10)
    labels = [c for c,_ in t]
    values = [v for _,v in t]
    plt.figure(figsize=(8,4))
    plt.bar(labels, values)
    plt.title(title)
    plt.xlabel("Codon"); plt.ylabel("Count")
    plt.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    plt.close()

def write_nutrition_outputs(outdir: Path, s2_aa3, flu_aa3):
    """
    Write a small nutrition dataset for the amino acids found in the top-3 lists.
    Produces:
      - amino_acid_foods.json
      - nutrition_summary.md
    """
    outdir.mkdir(parents=True, exist_ok=True)

    s2_top = [aa for aa,_ in s2_aa3]
    flu_top = [aa for aa,_ in flu_aa3]
    all_top = []
    for aa in s2_top + flu_top:
        if aa not in all_top:
            all_top.append(aa)

    aa_data = {}
    for aa in all_top:
        name = AA_NAMES.get(aa, aa)
        kb = NUTRITION_KB.get(aa, {
            "low_foods": ["Fruits", "Oils/fats"],
            "rich_foods": ["Animal proteins", "Dairy", "Legumes/soy"],
        })
        aa_data[aa] = {
            "name": name,
            "low_foods": kb["low_foods"],
            "rich_foods": kb["rich_foods"],
        }

    json_obj = {
        "sources": NUTRITION_SOURCES,
        "top_amino_acids": {
            "sars_cov_2": s2_top,
            "influenza": flu_top,
        },
        "amino_acids": aa_data,
    }
    (outdir / "amino_acid_foods.json").write_text(json.dumps(json_obj, indent=2))

    # Simple, human-readable summary
    lines = []
    lines.append("# Amino-acid Nutrition Summary\n")
    lines.append("**Top-3 (SARS-CoV-2):** " + ", ".join([f"{a} ({AA_NAMES.get(a,a)})" for a in s2_top]))
    lines.append("**Top-3 (Influenza):** " + ", ".join([f"{a} ({AA_NAMES.get(a,a)})" for a in flu_top]))
    lines.append("\n## Foods guidance (examples)\n")
    for aa in all_top:
        name = AA_NAMES.get(aa, aa)
        low_foods = ", ".join(NUTRITION_KB.get(aa, {}).get("low_foods", ["Fruits", "Oils"]))
        rich_foods = ", ".join(NUTRITION_KB.get(aa, {}).get("rich_foods", ["Animal proteins", "Dairy", "Legumes/soy"]))
        lines.append(f"### {aa} — {name}")
        lines.append(f"- Low in: {low_foods}")
        lines.append(f"- Rich in: {rich_foods}\n")
    lines.append("**Sources:** " + "; ".join(NUTRITION_SOURCES) + "\n")
    (outdir / "nutrition_summary.md").write_text("\n".join(lines))

def main():
    ap = argparse.ArgumentParser(description="Codon usage: SARS-CoV-2 vs Influenza")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--combined", help="Multi-FASTA containing both; use prefixes to split")
    src.add_argument("--pair", nargs=2, metavar=("COVID.fasta","FLU.fasta"),
                     help="Provide two separate FASTA files")

    ap.add_argument("--covid-prefix", default="covid|",
                    help="Prefix on headers in --combined belonging to SARS-CoV-2 (default: covid|)")
    ap.add_argument("--flu-prefix", default="influenza|",
                    help="Prefix on headers in --combined belonging to Influenza (default: influenza|)")
    ap.add_argument("--min-orf", type=int, default=90, help="Minimum ORF length (nt)")
    ap.add_argument("--outdir", default="labs_python/lab04/outputs",
                    help="Directory for charts and general outputs")
    ap.add_argument("--nutrition-outdir", default="labs_python/lab04/outputs/nutrition",
                    help="Directory to write nutrition data files")
    args = ap.parse_args()

    # Load sequences
    if args.combined:
        recs = list(read_fasta_records(args.combined))
        covid_recs = [(h,s) for h,s in recs if h.startswith(args.covid_prefix)]
        flu_recs   = [(h,s) for h,s in recs if h.startswith(args.flu_prefix)]
        if not covid_recs or not flu_recs:
            heads = [h for h,_ in recs]
            raise SystemExit(
                "Could not split combined FASTA.\n"
                f"Found headers:\n- " + "\n- ".join(heads) + "\n"
                f"Expected prefixes: '{args.covid_prefix}' and '{args.flu_prefix}'."
            )
        sars2_rna = cat_rna(covid_recs)
        flu_rna   = cat_rna(flu_recs)
    else:
        covid, flu = args.pair
        sars2_rna = cat_rna(read_fasta_records(covid))
        flu_rna   = cat_rna(read_fasta_records(flu))

    # ORFs -> codon counts
    sars2_orfs = find_orfs_rna(sars2_rna, min_len_nt=args.min_orf)
    flu_orfs   = find_orfs_rna(flu_rna,   min_len_nt=args.min_orf)

    sars2_cnt = count_codons(sars2_orfs)
    flu_cnt   = count_codons(flu_orfs)

    # Charts
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plot_top10(sars2_cnt, "Top 10 codons – SARS-CoV-2", str(outdir / "covid_top10_codons.png"))
    plot_top10(flu_cnt,   "Top 10 codons – Influenza A", str(outdir / "influenza_top10_codons.png"))

    # (c) Intersection of Top-10s
    sars2_top10 = [c for c,_ in sars2_cnt.most_common(10)]
    flu_top10   = [c for c,_ in flu_cnt.most_common(10)]
    common = [c for c in sars2_top10 if c in flu_top10]

    # (d) Top-3 amino acids for each genome
    s2_aa3  = aa_counts(sars2_cnt).most_common(3)
    flu_aa3 = aa_counts(flu_cnt).most_common(3)

    # Console output (no “prompt for AI”)
    print("\n=== Top 10 codons – SARS-CoV-2 ===")
    for c,v in sars2_cnt.most_common(10):
        print(f"{c}\t{v}\tAA:{GENETIC_CODE.get(c)}")

    print("\n=== Top 10 codons – Influenza A ===")
    for c,v in flu_cnt.most_common(10):
        print(f"{c}\t{v}\tAA:{GENETIC_CODE.get(c)}")

    print("\n=== Common codons in both Top-10 lists ===")
    print(", ".join(common) if common else "(none)")

    print("\n=== Top 3 amino acids – SARS-CoV-2 ===")
    print(", ".join([f"{aa}:{n}" for aa,n in s2_aa3]))
    print("=== Top 3 amino acids – Influenza A ===")
    print(", ".join([f"{aa}:{n}" for aa,n in flu_aa3]))

    # Write nutrition files for the AAs that actually appeared in the results
    write_nutrition_outputs(Path(args.nutrition_outdir), s2_aa3, flu_aa3)
    print(f"\nWrote nutrition data to: {args.nutrition_outdir}")

if __name__ == "__main__":
    main()
