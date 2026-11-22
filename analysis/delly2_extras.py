"""Two Delly2 comparisons quoted in the text (Discussion; Methods, Delly2; Additional file 1, Delly2 commands).

1. Per-sample Delly2: each genome's own PASS discovery calls (background 'none'), scored with the main-analysis
   filters and match rule, and the parents' calls longer than 1 Mb.
2. Germline workflow with no minimum site quality (-y 0, every Delly2 run of the paper) against Delly2's default
   minimum (-y 200): share of merged sites the minimum removes (sites.bcf of each merge, if kept:
   $CNVB_RESULTS/delly2_sites/{species}/{background}[_y200]/{sample}/sites.bcf), calls per individual that differ, filtered child deletions per
   background, and MIR under both. Needs bcftools on PATH for the site counts. Table delly2_extras.
    python analysis/delly2_extras.py"""
import subprocess

from cnvbench import calls as C, mir as M
from cnvbench.paths import CALLS, REPEATS, RESULTS, TABLES

BIG = 1_000_000
SITES = RESULTS / "delly2_sites"


def n_sites(path):
    out = subprocess.run(["bcftools", "view", "-H", str(path)], check=True, capture_output=True, text=True).stdout
    return out.count("\n")


def main():
    rows = []
    print("Per-sample Delly2 (discovery calls of each genome alone)")
    for sp in ("zea", "human"):
        tr = C.trios(sp)
        rep = M.load_repeats(sp, REPEATS[sp])
        calls = {s: C.load(sp, "Delly2", "none", s) for s in C.samples(sp)}
        head = M.apply_filters(calls, tr, rep)
        big = [sum(c.svtype == "DUP" and M.length(c) > BIG for c in head[p]) for t in tr for p in t.parents]
        for sv in M.SVTYPES:
            v = M.mir(head, tr, sv)
            rows.append((sp, "per_sample", "none", sv, "mir", f"{v:.1f}"))
            print(f"  {C.SPECIES[sp]} {sv} MIR {v:.1f}%")
        rows.append((sp, "per_sample", "none", "DUP", "parent_calls_over_1mb_mean", f"{sum(big) / len(big):.1f}"))
        print(f"  {C.SPECIES[sp]} filtered parental duplications > 1 Mb: {sum(big) / len(big):.1f} per parent "
              f"({min(big)}-{max(big)})")

    print("\nGermline workflow: -y 200 (Delly2 default) vs -y 0 (main)")
    for sp in ("zea", "human"):
        tr = C.trios(sp)
        rep = M.load_repeats(sp, REPEATS[sp])
        bgs = sorted(p.name for p in (CALLS / sp / "delly2").iterdir() if p.name != "none" and not p.name.endswith("_y200"))
        for bg in bgs:
            y200 = {s: C.load(sp, "Delly2", bg + "_y200", s) for s in C.samples(sp)}
            y0 = {s: C.load(sp, "Delly2", bg, s) for s in C.samples(sp)}
            diff = max(len(set(y200[s]) ^ set(y0[s])) for s in y200)
            h200, h0 = M.apply_filters(y200, tr, rep), M.apply_filters(y0, tr, rep)
            n_del = [sum(c.svtype == "DEL" for t in tr for c in h[t.child]) for h in (h200, h0)]
            mirs = {sv: (M.mir(h200, tr, sv), M.mir(h0, tr, sv)) for sv in M.SVTYPES}
            removed = []
            members = [p.name for p in (SITES / sp / bg).iterdir()] if (SITES / sp / bg).exists() else []
            for s in members:                    # one merge per trio member, or one ('all') for pooled-12
                a, b = SITES / sp / f"{bg}_y200" / s / "sites.bcf", SITES / sp / bg / s / "sites.bcf"
                if a.exists() and b.exists():
                    n200, n0 = n_sites(a), n_sites(b)
                    removed.append(100 * (n0 - n200) / n0)
            rm = f"{min(removed):.0f}-{max(removed):.0f} (mean {sum(removed) / len(removed):.0f})" if removed else "n.a."
            rows.append((sp, "y200_vs_y0", bg, "", "max_calls_differing_per_individual", diff))
            rows.append((sp, "y200_vs_y0", bg, "DEL", "filtered_child_calls_y200_y0", f"{n_del[0]}/{n_del[1]}"))
            rows.append((sp, "y200_vs_y0", bg, "", "merged_sites_removed_pct", rm))
            for sv, (a, b) in mirs.items():
                rows.append((sp, "y200_vs_y0", bg, sv, "mir_y200_y0", f"{a:.1f}/{b:.1f}"))
            print(f"  {C.SPECIES[sp]} {bg:14s} sites removed by -y 200: {rm}% (per member); calls differing per individual "
                  f"<= {diff}; filtered child deletions {n_del[0]} (y0 {n_del[1]}); "
                  + "; ".join(f"{sv} MIR {a:.1f} vs {b:.1f}" for sv, (a, b) in mirs.items()))

    TABLES.mkdir(parents=True, exist_ok=True)
    with open(TABLES / "delly2_extras.tsv", "w") as fh:
        fh.write("species\tcomparison\tbackground\tsvtype\tmeasure\tvalue\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    print("wrote", TABLES / "delly2_extras.tsv")


if __name__ == "__main__":
    main()
