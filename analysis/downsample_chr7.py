"""Chromosome 7 depth series (Additional file 1, Section S12): the 12 trio moths at full depth, 10x and 5x
(workflow/downsample), each caller with its headline background at the same depth. Copy number at the Tryp77
amplicon against ddPCR, and chromosome 7 MIR of the four children under the headline match rule and filters.
Tables: downsample_chr7_summary (Table S21), downsample_chr7_cn, downsample_chr7_depth. Figure: downsample_chr7
(Fig. S19).
    python analysis/downsample_chr7.py
Each level is its own results tree (results/downsample/LEVEL), and cnvbench reads CNVB_RESULTS when it is
imported, so every level is scored in a child process (downsample_chr7.py level LEVEL) that caches
results/downsample/LEVEL/tables/downsample_chr7.json."""
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

LEVELS = ["full", "10x", "5x"]
CALLERS = ["Delly2", "smoove", "GATK-gCNV", "CNVcaller", "CNVkit"]
CN_CALLERS = [c for c in CALLERS if c != "smoove"]      # smoove reports no copy number
CN_DESIGN = {"Delly2": "crosstrio9", "smoove": "crosstrio9", "GATK-gCNV": "crosstrio9", "CNVcaller": "crosstrio9",
             "CNVkit": "crosstrio3"}        # Delly2: population delly cnv workflow (delly2_pop), as at full depth


def level():
    """Copy number and chromosome 7 MIR at the level whose tree CNVB_RESULTS points to."""
    from cnvbench import calls as C, locus as L, mir as M
    from cnvbench.paths import REPEATS, TABLES, call_dir
    trios, rep = C.trios("zea"), M.load_repeats("zea", REPEATS["zea"])
    res = {"cn": {}, "mir": {}}
    for m in CALLERS:
        if m in CN_CALLERS:
            res["cn"][m] = {s: L.caller_cn(m, CN_DESIGN[m], s) for s in C.samples("zea")}
        raw = {s: [c for c in C.load("zea", m, C.HEADLINE[m], s) if c.chrom == L.CHROM] for s in C.samples("zea")}
        filt = M.apply_filters(raw, trios, rep)
        res["mir"][m] = {}
        for sv in M.SVTYPES:
            pt = M.per_trio(filt, trios, sv)
            ptc = M.per_trio(filt, trios, sv, parent_max=1_000_000)       # parental calls > 1 Mb may not match
            res["mir"][m][sv] = {"n": sum(n for _, n in pt), "mir": M.pct(sum(k for k, _ in pt), sum(n for _, n in pt)),
                                 "mir_le1mb": M.pct(sum(k for k, _ in ptc), sum(n for _, n in ptc))}
    # quality of each moth's delly cnv call over the ddPCR amplicon (no breakpoint refinement, as run here);
    # delly merge -y filters on this quality in the monitoring workflow
    res["delly_qual"] = {}
    for s in C.samples("zea"):
        out = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%INFO/END\t%QUAL\n",
                              str(call_dir("zea", "delly2_cnv", "none", s) / f"{s}.bcf")],
                             capture_output=True, text=True, check=True).stdout.split("\n")
        q = [float(x[3]) for x in (l.split("\t") for l in out if l)
             if C.norm_chrom("zea", x[0]) == L.CHROM and int(x[1]) <= L.AMPLICON[1] and int(x[2]) >= L.AMPLICON[0]]
        res["delly_qual"][s] = max(q, default=None)
    TABLES.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(TABLES / "downsample_chr7.json", "w"), indent=1)


def summary():
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import stats

    from cnvbench import locus as L
    from cnvbench.paths import RESULTS, TABLES
    from cnvbench.style import COLORS, MM, label, save, tint
    R = {}
    for lev in LEVELS:
        cache = RESULTS / "downsample" / lev / "tables" / "downsample_chr7.json"
        if not cache.exists():
            subprocess.run([sys.executable, __file__, "level", lev], check=True,
                           env={**os.environ, "CNVB_RESULTS": str(RESULTS / "downsample" / lev)})
        R[lev] = json.load(open(cache))
    dd = L.ddpcr()
    moths = list(dd)
    TABLES.mkdir(parents=True, exist_ok=True)
    f = lambda v: "" if v is None else f"{v:.3f}" if isinstance(v, float) else str(v)
    with open(TABLES / "downsample_chr7_summary.tsv", "w") as fh, open(TABLES / "downsample_chr7_cn.tsv", "w") as fc:
        fh.write("caller\tdepth\tr\trho\tcarriers_called\tfalse_noncarriers\tchild_dup\tchild_del\tdup_mir\tdel_mir\t"
                 "dup_mir_le1mb\tdel_mir_le1mb\n")
        fc.write("depth\tcaller\tsample\tddpcr\tcalled\n")
        for m in CALLERS:
            for lev, res in R.items():
                cn, mr = res["cn"].get(m), res["mir"][m]
                car = [s for s in moths if dd[s] >= 2.5]
                non = [s for s in moths if dd[s] < 2.5]
                if cn:
                    x, y = [dd[s] for s in moths], [cn[s] for s in moths]
                    cnv = [stats.pearsonr(x, y)[0], stats.spearmanr(x, y)[0], f"{sum(cn[s] >= 2.5 for s in car)}/{len(car)}",
                           f"{sum(cn[s] >= 2.5 or cn[s] <= 1.5 for s in non)}/{len(non)}"]
                else:
                    cnv = [None] * 4
                row = [m, lev, *cnv,
                       mr["DUP"]["n"], mr["DEL"]["n"], mr["DUP"]["mir"], mr["DEL"]["mir"],
                       mr["DUP"]["mir_le1mb"], mr["DEL"]["mir_le1mb"]]
                fh.write("\t".join(f(v) for v in row) + "\n")
                if cn:
                    fc.writelines(f"{lev}\t{m}\t{s}\t{dd[s]:.2f}\t{cn[s]:.2f}\n" for s in moths)
    with open(TABLES / "downsample_chr7_depth.tsv", "w") as fh:
        fh.write("sample\tlevel\tfraction\tdepth_bt2\tdepth_bwa\tfull_depth_bt2\n")
        for s in moths:
            for r in list(csv.reader(open(RESULTS / "downsample" / "depth" / f"{s}.tsv"), delimiter="\t"))[1:]:
                fh.write(s + "\t" + "\t".join(r) + "\n")
    print(open(TABLES / "downsample_chr7_summary.tsv").read())
    with open(TABLES / "downsample_chr7_delly_qual.tsv", "w") as fh:
        fh.write("sample\tddpcr\t" + "\t".join(LEVELS) + "\n")
        for s in sorted(moths, key=lambda s: -dd[s]):
            q = [R[lev]["delly_qual"][s] for lev in LEVELS]
            fh.write(f"{s}\t{dd[s]:.2f}\t" + "\t".join("" if v is None else f"{v:.0f}" for v in q) + "\n")
    print("Delly2 cnv quality at the amplicon (full, 10x, 5x):\n" + open(TABLES / "downsample_chr7_delly_qual.tsv").read())

    # Fig. S19: a copy number against ddPCR per caller and level; b chromosome 7 MIR against depth
    shade = {"full": 1.0, "10x": 0.6, "5x": 0.3}
    mk = {"full": "o", "10x": "s", "5x": "D"}
    name = {"full": "Full depth", "10x": "10×", "5x": "5×"}
    fig = plt.figure(figsize=(180 * MM, 110 * MM))
    gs = fig.add_gridspec(2, 5, height_ratios=[1, 1.05], hspace=0.55, wspace=0.38, left=0.06, right=0.99, top=0.93,
                          bottom=0.1)
    for j, m in enumerate(CN_CALLERS):
        ax = fig.add_subplot(gs[0, j])
        ax.plot([1.5, 10.5], [1.5, 10.5], color="#BBBBBB", lw=0.6, ls=":", zorder=1)
        for lev, res in R.items():
            c = COLORS[m]
            ax.scatter([dd[s] for s in moths], [res["cn"][m][s] for s in moths], s=12, marker=mk[lev], lw=0.6,
                       facecolor=tint(c, shade[lev] * 0.55), edgecolor=tint(c, shade[lev]),
                       zorder=3 if lev == "full" else 2, label=name[lev])
        ax.set(xlim=(1.3, 9.2), ylim=(-0.3, 10.6), xticks=[2, 4, 6, 8], yticks=[0, 2, 4, 6, 8, 10])
        ax.set_title(m, loc="left")
        ax.grid(True, zorder=0); ax.set_axisbelow(True)
        ax.set_xlabel("ddPCR copy number")
        if j == 0:
            ax.set_ylabel("Called copy number")
            label(ax, "a", x=-0.38)
    key = fig.add_subplot(gs[0, 4])
    key.axis("off")
    key.legend(*ax.get_legend_handles_labels(), frameon=False, loc="center left", handletextpad=0.3)
    x = np.arange(len(R))
    for k, sv in enumerate(("DUP", "DEL")):
        ax = fig.add_subplot(gs[1, 2 * k:2 * k + 2])
        for m in CALLERS[:-1]:              # CNVkit: 3-5 child calls on chromosome 7, too few for a rate
            y = [R[lev]["mir"][m][sv]["mir"] for lev in R]
            ax.plot(x, [np.nan if v is None else v for v in y], color=COLORS[m], lw=1.3, marker="o", ms=3.5, label=m,
                    ls="--" if m in ("Delly2", "smoove") else "-")
        y = [R[lev]["mir"]["smoove"][sv]["mir_le1mb"] for lev in R]
        ax.plot(x, [np.nan if v is None else v for v in y], color=COLORS["smoove"], lw=0.9, ls=":", marker="o", ms=3,
                mfc="white", label="smoove, parental calls\n> 1 Mb may not match")
        ax.set_xticks(x, [{"full": "Full", "10x": "10×", "5x": "5×"}[lev] for lev in R])
        ax.set(xlabel="Chromosome 7 depth", ylabel="MIR (%)" if k == 0 else None, xlim=(-0.3, len(R) - 0.7))
        ax.set_ylim(bottom=0)
        ax.set_title(f"Chromosome 7 {'duplications' if sv == 'DUP' else 'deletions'}", loc="left")
        ax.grid(axis="y", zorder=0); ax.set_axisbelow(True)
        if k == 0:
            label(ax, "b", x=-0.17)
    fig.legend(*ax.get_legend_handles_labels(), loc="center left", bbox_to_anchor=(0.815, 0.27), frameon=False)
    save(fig, "downsample_chr7", max_mm=None)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    level() if sys.argv[1:2] == ["level"] else summary()
