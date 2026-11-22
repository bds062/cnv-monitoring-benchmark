"""Copy number at Tryp77 against ddPCR: caller copy number of the 12 trio moths under every background,
correlations (Additional file 1, Tables S17 and S24), Fig. 3 (fig_ddpcr_panels) and Fig. S25 (locus_depth_calls).
    python analysis/ddpcr.py"""
import json
import random

import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle
from scipy import stats

from cnvbench import calls as C, figstyle, locus as L, mir as M, style
from cnvbench.paths import REPEATS, TABLES, bam
from cnvbench.style import COLORS, GREY, INK, MM, MUTED, ORDER

SP = "zea"
TRIOS = C.trios(SP)
ALL = C.samples(SP)
ROLE = C.roles(SP)
DD = L.ddpcr()
# copy number at Tryp77 under every background; Delly2's is its population read-depth workflow (delly cnv refined by
# breakpoints, delly merge -y 0, delly cnv -u, delly classify), the Delly2 workflow of the monitoring analyses
CN_HEADLINE = {**C.HEADLINE, "Delly2": "crosstrio9"}
BACKGROUNDS = {"Delly2": ["pooled12", "susceptible13", "neutral18", "neutral9", "crosstrio9"],
               "CNVcaller": ["pooled12", "susceptible13", "neutral18", "neutral9", "crosstrio9"],
               "GATK-gCNV": ["pooled12", "susceptible13", "neutral18", "neutral9", "crosstrio9"],
               "CNVkit": ["pooled12", "susceptible13", "neutral18", "crosstrio3"]}
CN_CALLERS = [m for m in ORDER if m != "smoove"]   # smoove reports no copy number, so it is not compared with ddPCR
HIGH = [s for s in ALL if DD[s] > 6]             # the three resistant parents above six copies
MARK = {"resistant": "^", "susceptible": "v", "child": "o"}
ROLE_LABEL = {"resistant": "Resistant parent", "susceptible": "Susceptible parent", "child": "Child"}


def corr(vals, samples=ALL):
    x, y = [DD[s] for s in samples], [vals[s] for s in samples]
    if len(set(y)) < 2:
        return np.nan, np.nan
    return stats.pearsonr(x, y)[0], stats.spearmanr(x, y)[0]


def boot_r(vals, rng, reps=4000):
    """95% interval of Pearson r over resamples of the 12 moths."""
    x, y = np.array([DD[s] for s in ALL]), np.array([vals[s] for s in ALL])
    bs = []
    for _ in range(reps):
        i = [rng.randrange(len(x)) for _ in x]
        if len(set(x[i])) > 2 and len(set(y[i])) > 1:
            bs.append(stats.pearsonr(x[i], y[i])[0])
    return np.percentile(bs, [2.5, 97.5]).tolist()


def tables():
    cn = {m: {s: L.caller_cn(m, CN_HEADLINE[m], s) for s in ALL} for m in CN_CALLERS}
    with open(TABLES / "tryp77_cn.tsv", "w") as fh:            # Additional file 1, Table S17
        fh.write("trio\trole\tsample\tddpcr\t" + "\t".join(CN_CALLERS) + "\n")
        for t in TRIOS:
            for s in (*t.parents, t.child):
                fh.write(f"{t.name}\t{ROLE[s]}\t{s}\t{DD[s]:.2f}\t" + "\t".join(f"{cn[m][s]:.2f}" for m in CN_CALLERS) + "\n")
    rng, rows = random.Random(7), []
    for m in ["Delly2", "GATK-gCNV", "CNVcaller", "CNVkit"]:
        r, rho = corr(cn[m])
        lo, hi = boot_r(cn[m], rng)
        r_low = corr(cn[m], [s for s in ALL if DD[s] <= 6])[0]
        rows.append((m, r, rho, lo, hi, r_low))
    with open(TABLES / "tryp77_ddpcr_correlation.tsv", "w") as fh:
        fh.write("caller\tpearson_r\tspearman_rho\tr_boot_lo\tr_boot_hi\tr_without_ddpcr_gt6\n")
        fh.writelines("\t".join([m] + [f"{v:.3f}" for v in vals]) + "\n" for m, *vals in rows)
    with open(TABLES / "tryp77_ddpcr_by_background.tsv", "w") as fh:     # Additional file 1, Table S24
        fh.write("caller\tbackground\theadline\tpearson_r\tspearman_rho\tcn_ddpcr_gt6\n")
        for m, bgs in BACKGROUNDS.items():
            for bg in bgs:
                v = cn[m] if bg == CN_HEADLINE[m] else {s: L.caller_cn(m, bg, s) for s in ALL}
                r, rho = corr(v)
                fh.write(f"{m}\t{bg}\t{int(bg == CN_HEADLINE[m])}\t{r:.3f}\t{rho:.3f}\t"
                         + ", ".join(f"{v[s]:.2f}" for s in HIGH) + "\n")
    return cn


def fig_ddpcr(cn):
    """Fig. 3 (print style, 170 mm): one panel per caller, caller copy number against ddPCR."""
    F = figstyle
    with F.paper():
        W, H = F.FULL, 54.0
        fig = plt.figure(figsize=(W * F.MM, H * F.MM))
        L_, R_, GAP, B_, T_ = 13.0, 1.0, 5.0, 15.5, 5.5                  # margins and gap (mm)
        w = (W - L_ - R_ - 3 * GAP) / 4
        h = min(w, H - B_ - T_)
        lim = (0, 9.5)
        for k, m in enumerate(CN_CALLERS):
            ax = fig.add_axes([(L_ + k * (w + GAP)) / W, B_ / H, w / W, h / H])
            ax.plot(lim, lim, ls=":", color="black", lw=F.LW, zorder=1)
            for role, mk in MARK.items():
                ss = [s for s in ALL if ROLE[s] == role]
                ax.scatter([DD[s] for s in ss], [cn[m][s] for s in ss], marker=mk, s=30, zorder=3,
                           facecolor=F.COLORS[m], edgecolor="black", lw=F.THIN, alpha=0.9)
            r, rho = corr(cn[m])
            ax.text(0.05, 0.95, f"$r$ = {r:.2f}\n$\\rho$ = {rho:.2f}", transform=ax.transAxes, va="top",
                    fontsize=F.SMALL, linespacing=1.3)
            ax.set(xlim=lim, ylim=lim, xticks=[0, 2, 4, 6, 8], yticks=[0, 2, 4, 6, 8])
            ax.grid(True, axis="both")
            ax.tick_params(labelleft=k == 0)
            ax.set_title(m, loc="left")
            if k == 0:
                ax.set_ylabel("Caller copy number")
            F.label(ax, "abcd"[k], x=-11.5 / w if k == 0 else -4.6 / w, y=1.0)
        fig.text((L_ + (W - L_ - R_) / 2) / W, (B_ - 7.2) / H, "ddPCR copy number", ha="center", va="center",
                 fontsize=F.FS)
        hs = [plt.Line2D([], [], ls="", marker=mk, ms=5.5, mfc=F.GREY, mec="black", mew=F.THIN)
              for mk in MARK.values()] + [plt.Line2D([], [], ls=":", color="black", lw=F.LW)]
        fig.legend(hs, [ROLE_LABEL[r] for r in MARK] + ["$y$ = $x$"], loc="lower center", ncol=4,
                   bbox_to_anchor=((L_ + (W - L_ - R_) / 2) / W, 0.3 / H))
        F.save(fig, "fig_ddpcr_panels")


def locus_depth():
    """Read depth over the amplicon +- 100 kb in 1 kb windows, as 2 x depth / the moth's median depth over
    chromosome 7 in 50 kb windows (BWA-MEM alignments); cached in tryp77_depth.json."""
    cache = TABLES / "tryp77_depth.json"
    if cache.exists():
        return json.load(open(cache))
    w0, w1 = L.AMPLICON[0] - 100_000, L.AMPLICON[1] + 100_000
    dep = {"windows": [(a, a + 1000) for a in range(w0 // 1000 * 1000, w1, 1000)]}
    for s in ALL:
        path = bam(SP, s, "bam_bwa")
        contig = L.bam_contig(path, SP, L.CHROM)
        n = pysam.AlignmentFile(str(path)).get_reference_length(contig)
        med = float(np.median(L.bam_depth(path, contig, [(a, a + 50_000) for a in range(0, n - 50_000, 50_000)])))
        dep[s] = [2 * d / med for d in L.bam_depth(path, contig, dep["windows"])]
    json.dump(dep, open(cache, "w"))
    return dep


def fig_locus():
    dep = locus_depth()
    w0, w1 = L.AMPLICON[0] - 100_000, L.AMPLICON[1] + 100_000
    xs = np.array([(a + b) / 2e6 for a, b in dep["windows"]])
    smooth = lambda v: np.convolve(v, np.ones(3) / 3, mode="same")
    susc = np.mean([dep[s] for s in ALL if ROLE[s] == "susceptible"], axis=0)
    raw = C.headline(SP)
    rep = M.load_repeats(SP, REPEATS[SP])
    kept = {m: M.apply_filters(raw[m], TRIOS, rep) for m in ORDER}
    lane = lambda k: -k - (0.5 if k >= 2 else 0)            # gap between breakpoint and read-depth callers
    fig = plt.figure(figsize=(170 * MM, 160 * MM))
    outer = fig.add_gridspec(3, 4, hspace=0.14, wspace=0.12, left=0.11, right=0.99, top=0.965, bottom=0.105)
    for col, t in enumerate(TRIOS):
        for row, s in enumerate((*t.parents, t.child)):
            sub = outer[row, col].subgridspec(2, 1, height_ratios=[1.0, 0.85], hspace=0.06)
            ax, al = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
            for a_ in (ax, al):
                for _, a, b in L.GENES:
                    a_.axvspan(a / 1e6, b / 1e6, color="#EDEDED", lw=0, zorder=0)
                a_.axvline(sum(L.AMPLICON) / 2e6, color=GREY, lw=0.6, ls=":", zorder=1)
                a_.set_xlim(w0 / 1e6, w1 / 1e6)
            ax.plot(xs, smooth(susc), color="#B5B5B5", lw=0.8)
            ax.plot(xs, np.clip(smooth(np.array(dep[s])), 0, 10), color=INK, lw=0.8)
            ax.axhline(2, color="#DDDDDD", lw=0.5, zorder=0)
            ax.set(ylim=(0, 10), yticks=[0, 2, 4, 6, 8, 10])
            ax.tick_params(labelbottom=False, labelleft=col == 0)
            if col == 0:
                ax.set_ylabel("Depth CN")
                ax.text(-0.42, 0.5, ROLE_LABEL[ROLE[s]], transform=ax.transAxes, rotation=90, ha="center",
                        va="center", fontsize=7.5, fontweight="bold")
            if row == 0:
                ax.set_title(t.name, loc="left", fontweight="bold")
            ax.text(0.98, 0.95, f"ddPCR {DD[s]:.1f}", transform=ax.transAxes, ha="right", va="top", fontsize=6.3,
                    fontweight="bold" if DD[s] > 2.5 else None, color=INK if DD[s] > 2.5 else MUTED)
            for k, m in enumerate(ORDER):
                ok = set(kept[m][s])
                cs = [c for c in raw[m][s] if c.chrom == L.CHROM and c.svtype == "DUP" and c.end > w0 and c.start < w1]
                for z, c in enumerate(sorted(cs, key=lambda c: c.start - c.end)):     # longest first, nested on top
                    a, b = max(c.start, w0) / 1e6, min(c.end, w1) / 1e6
                    h = 0.24 if c.end - c.start > 1_000_000 else 0.7
                    al.add_patch(Rectangle((a, lane(k) - h / 2), max(b - a, 0.0008), h, zorder=2 + 0.01 * z,
                                           facecolor=to_rgb(COLORS[m]) + (0.75,) if c in ok else "white",
                                           edgecolor=COLORS[m], lw=0.6 if c in ok else 0.8))
            al.axhline(lane(1) - 0.75, color="#BBBBBB", lw=0.6, ls=":")
            al.set_ylim(lane(4) - 0.6, 0.6)
            al.set_yticks([lane(k) for k in range(5)], ORDER if col == 0 else [""] * 5, fontsize=6.3)
            al.tick_params(axis="y", length=0)
            al.spines["left"].set_visible(False)
            al.set_xticks([8.20, 8.25, 8.30])
            al.tick_params(labelbottom=row == 2)
    grey = to_rgb("#666666") + (0.75,)
    h = [plt.Line2D([], [], color=INK, lw=0.8, label="This moth"),
         plt.Line2D([], [], color="#B5B5B5", lw=0.8, label="Susceptible parents (mean)"),
         Rectangle((0, 0), 1, 1, facecolor=grey, edgecolor="#666666", lw=0.6, label="Duplication call"),
         Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="#666666", lw=0.8, label="Removed by filters"),
         plt.Line2D([], [], color=grey, lw=2.2, solid_capstyle="butt", label="Duplication > 1 Mb"),
         Rectangle((0, 0), 1, 1, facecolor="#EDEDED", lw=0, label="Tryp079 / 078 / 077"),
         plt.Line2D([], [], color=GREY, lw=0.6, ls=":", label="ddPCR amplicon")]
    fig.text(0.55, 0.056, "Chromosome 7 (Mb)", ha="center", va="center", fontsize=7)
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.55, -0.022),
               handlelength=1.6, columnspacing=1.2)
    style.save(fig, "locus_depth_calls")


if __name__ == "__main__":
    TABLES.mkdir(parents=True, exist_ok=True)
    fig_ddpcr(tables())
    fig_locus()
