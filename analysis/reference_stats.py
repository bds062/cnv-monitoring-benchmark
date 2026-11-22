"""Reference facts quoted in Methods: the share of each reference masked by the repeat track used in the repeat
filter (H. zea: RefSeq soft-masking, chromosomes 1-31; human: UCSC RepeatMasker, chromosomes 1-22), and the
lengths of the lifted OGS2 gene models that set the MIR size bins (data/ogs2_gene_lengths.tsv). Table reference_stats.
    python analysis/reference_stats.py"""
import numpy as np

from cnvbench import calls as C, mir as M
from cnvbench.paths import REF, REPEATS, REPO, TABLES

CHROMS = {"zea": {str(i) for i in range(1, 32)}, "human": {str(i) for i in range(1, 23)}}


def main():
    rows = []
    for sp in ("zea", "human"):
        rep = M.load_repeats(sp, REPEATS[sp])
        size = {C.norm_chrom(sp, l.split("\t")[0]): int(l.split("\t")[1]) for l in open(f"{REF[sp]}.fai")}
        keep = [c for c in CHROMS[sp] if c in size]
        masked = sum(rep[c][2][-1] for c in keep if c in rep)
        total = sum(size[c] for c in keep)
        rows.append((sp, "repeat_masked_pct", f"{100 * masked / total:.1f}"))
        print(f"{C.SPECIES[sp]}: {100 * masked / total:.1f}% of chromosomes {min(map(int, keep))}-{max(map(int, keep))} masked")

    lengths = np.array([int(r.split("\t")[-1]) for r in open(REPO / "data" / "ogs2_gene_lengths.tsv")
                        if not r.startswith("#") and r.split("\t")[-1].strip().isdigit()])
    shares = [100 * np.mean((lengths >= lo) & (lengths < hi)) for lo, hi, _ in M.BINS]
    stats = {"genes": len(lengths), "min_bp": lengths.min(), "max_bp": lengths.max(), "median_bp": np.median(lengths),
             "mean_bp": lengths.mean(), "sd_bp": lengths.std(ddof=1)}
    for k, v in stats.items():
        rows.append(("zea", f"gene_{k}", f"{v:.0f}"))
    for (_, _, lab), v in zip(M.BINS, shares):
        rows.append(("zea", f"genes_in_{lab}_pct", f"{v:.1f}"))
    print(f"OGS2 genes placed: {len(lengths):,}; length {lengths.min()} bp - {lengths.max() / 1e3:.0f} kb, "
          f"median {stats['median_bp']:,.0f}, mean {stats['mean_bp']:,.0f}, SD {stats['sd_bp']:,.0f} bp; "
          f"bins {', '.join(f'{lab} {v:.1f}%' for (_, _, lab), v in zip(M.BINS, shares))} (genes < 100 bp outside the bins)")
    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "reference_stats.tsv", "w") as fh:
        fh.write("species\tmeasure\tvalue\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)
    print("wrote", TABLES / "reference_stats.tsv")


if __name__ == "__main__":
    main()
