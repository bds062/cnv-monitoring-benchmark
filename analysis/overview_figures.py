"""Study-design schematic and caller scorecard, drawn as pictures rather than text.
  study_design      one column per evaluation criterion (a trio Mendelian consistency, b trio copy number,
                    c copy number in field monitoring): data -> the callers by kind of evidence -> how each is scored (no results)
  caller_scorecard  callers x criteria; each cell is a bar on that column's own scale (value underneath, best in bold),
                    monitoring outcomes as moth dots and false-call marks
    python overview_figures.py
Reads the tables of trio_benchmark.py, ddpcr.py, human_truthset.py, monitoring_panels.py and monitoring_scan.py;
run those first."""
import csv
import math
from pathlib import Path

import matplotlib.patches as mp
import matplotlib.pyplot as plt
import numpy as np

from cnvbench import calls as C, figstyle, locus as L
from cnvbench.paths import TABLES
import runtime
from cnvbench.style import COLORS, GREY as DARK, INK, MM, MUTED, ORDER, save, tint

PANEL = "#F3F3F3"                           # column backgrounds (no caller hue)
SYM = dict(fontfamily="DejaVu Sans")        # check, cross and arrow glyphs
SUB = "#5A5A5A"                             # explanatory labels: darker than MUTED
MON = ["GATK-gCNV", "CNVcaller", "Delly2"]  # callers run on the time series
MON_ROW = {"GATK-gCNV": "GATK-gCNV", "CNVcaller": "CNVcaller", "Delly2": "Delly2 merge (-y 0)"}


def tsv(name):
    return list(csv.DictReader(open(TABLES / f"{name}.tsv"), delimiter="\t"))


DD = L.ddpcr()
TRIOS = {t.name: tuple(DD[s] for s in (*t.parents, t.child)) for t in C.trios("zea")}   # resistant, susceptible, child
PC = {(r["panel"], r["caller"]): r for r in tsv("monitoring_panel_calls")}
SCAN = {r["caller"]: r for r in tsv("monitoring_scan")}
RARE = {m: int(PC["2002+2012", MON_ROW[m]]["carriers_called"]) for m in MON}       # the 1 carrier in 24
COMMON = {m: int(PC["2017", MON_ROW[m]]["carriers_called"]) for m in MON}          # carriers called of 8 in 11
FLAGGED = {m: int(SCAN[MON_ROW[m]]["flagged"]) for m in MON}                       # chr 7 sites rising 2002 -> 2017


# ── shared drawing helpers (mm coordinates, y up) ─────────────────────────────
def canvas(w, h):
    fig = plt.figure(figsize=(w * MM, h * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, w); ax.set_ylim(0, h); ax.axis("off")
    return fig, ax


def txt(ax, x, y, s, size=6.2, ha="left", va="center", color=INK, **kw):
    ax.text(x, y, s, fontsize=size, ha=ha, va=va, color=color, linespacing=1.15, **kw)


def box(ax, x, y, w, h, fc=PANEL, ec="none", r=1.6, lw=0.6):
    ax.add_patch(mp.FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw))


def arrow(ax, x0, y0, x1, y1, color=MUTED, lw=0.8):
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=7,
                                                         shrinkA=0, shrinkB=0))


def chip(ax, x, y, m, size=6.2):
    ax.add_patch(mp.Circle((x, y), 0.9, fc=COLORS[m], ec="none"))
    txt(ax, x + 1.7, y, m, size=size)


def dot(ax, x, y, filled, r=0.75, color=DARK):
    ax.add_patch(mp.Circle((x, y), r, fc=color if filled else "white", ec=color, lw=0.5, zorder=3))


def stack(ax, x, y, cn, s=0.72, gap=0.2, w=1.25):
    """Tryp77 copies (ddPCR, rounded) as a column of blocks rising from y."""
    for i in range(int(round(cn))):
        ax.add_patch(mp.Rectangle((x, y + i * (s + gap)), w, s, fc=DARK, ec="none"))


def ind(ax, x, y, sex, r=1.6, fc="white", ec=DARK, ls="-"):
    """Pedigree symbol: square for a male (ZZ), circle for a female (ZW)."""
    if sex == "M":
        ax.add_patch(mp.Rectangle((x - r * 0.9, y - r * 0.9), r * 1.8, r * 1.8, fc=fc, ec=ec, lw=0.6, ls=ls, zorder=3))
    else:
        ax.add_patch(mp.Circle((x, y), r, fc=fc, ec=ec, lw=0.6, ls=ls, zorder=3))


def mini_trio(ax, x, y, cns, sexes, dx=11, dy=8, r=1.6):
    (cr, cs, cc) = cns
    ax.plot([x, x + dx], [y, y], color=DARK, lw=0.6); ax.plot([x + dx / 2] * 2, [y, y - dy], color=DARK, lw=0.6)
    for (xx, yy), sx in zip(((x, y), (x + dx, y), (x + dx / 2, y - dy)), sexes):
        ind(ax, xx, yy, sx, r)
    stack(ax, x - r - 1.9, y - r, cr); stack(ax, x + dx + r + 0.6, y - r, cs); stack(ax, x + dx / 2 + r + 0.6, y - dy - r, cc)


# ── 1. study design (print style, 170 x 112 mm) ───────────────────────────────
# Three columns, one per evaluation criterion; each runs data (top) -> the callers (one strip shared by the
# three: breakpoint | read depth) -> evaluation (bottom). Schematic only: no results are drawn. The dots beside
# each evaluation title are the callers it covers.
SEX = {"Trio1": "FMF", "Trio2": "MFF", "Trio3": "FMF", "Trio4": "MFM"}    # resistant, susceptible, child (Z depth)
CN_CALLERS = ["Delly2", "GATK-gCNV", "CNVcaller", "CNVkit"]                # report copy number


def study_design():
    S = figstyle
    INKB, MID, FAINT = "black", "#595959", "#bdbdbd"
    T1, T2, T3 = S.FS, 8.5, S.SMALL                         # column heads, box titles, labels (pt)
    W, H = S.FULL, 112.0
    CW, CX = 52.5, (7.0, 61.25, 115.5)                      # column width, column left edges (mm)
    YD, YC, YE = (68.0, 103.0), (43.0, 62.0), (0.5, 37.0)   # data, callers, evaluation rows (bottom, top)
    with S.paper():
        fig, ax = canvas(W, H)

        def t(x, y, s, size=T3, **kw):
            kw.setdefault("ha", "left"); kw.setdefault("va", "center")
            ax.text(x, y, s, fontsize=size, color=kw.pop("color", INKB), linespacing=1.1, **kw)

        def frame(x, y, w, h, fc="white", ec=INKB, lw=S.LW, z=0):
            ax.add_patch(mp.Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=z))

        def arr(x0, y0, x1, y1):
            ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=INKB, lw=S.LW,
                                                                 mutation_scale=8, shrinkA=0, shrinkB=0))

        def person(x, y, sex, r=1.4, fc="white", ec=INKB, ls="-"):
            """Pedigree symbol: square for a male (ZZ), circle for a female (ZW)."""
            if sex == "M":
                ax.add_patch(mp.Rectangle((x - r * 0.9, y - r * 0.9), r * 1.8, r * 1.8, fc=fc, ec=ec, lw=S.THIN + 0.1,
                                          ls=ls, zorder=3))
            else:
                ax.add_patch(mp.Circle((x, y), r, fc=fc, ec=ec, lw=S.THIN + 0.1, ls=ls, zorder=3))

        def copies(x, y, cn, s=0.72, gap=0.2, w=1.3):
            """Tryp77 copies (ddPCR, rounded) as a column of blocks rising from y."""
            for i in range(int(round(cn))):
                ax.add_patch(mp.Rectangle((x, y + i * (s + gap)), w, s, fc=MID, ec="none"))

        def trio(x, y, cns, sexes, dx=11.0, dy=6.0, r=1.4):
            ax.plot([x, x + dx], [y, y], color=INKB, lw=S.THIN + 0.1)
            ax.plot([x + dx / 2] * 2, [y, y - dy], color=INKB, lw=S.THIN + 0.1)
            for (xx, yy), sx in zip(((x, y), (x + dx, y), (x + dx / 2, y - dy)), sexes):
                person(xx, yy, sx, r)
            copies(x - r - 1.9, y - r, cns[0]); copies(x + dx + r + 0.6, y - r, cns[1])
            copies(x + dx / 2 + r + 0.6, y - dy - r, cns[2])

        def callers(xr, y, ms):
            """Caller dots ending at xr: the callers an evaluation covers."""
            for k, m in enumerate(ms[::-1]):
                ax.add_patch(mp.Circle((xr - k * 2.7, y), 1.0, fc=S.COLORS[m], ec=INKB, lw=0.4, zorder=3))

        def tick(x, y, ok, s=1.3):
            """Check (ok) or cross mark drawn as lines."""
            if ok:
                ax.plot([x - s, x - 0.3 * s, x + s], [y, y - 0.8 * s, y + 0.9 * s], color=INKB, lw=1.1,
                        solid_capstyle="round", solid_joinstyle="round")
            else:
                for d in (1, -1):
                    ax.plot([x - 0.8 * s, x + 0.8 * s], [y - d * 0.8 * s, y + d * 0.8 * s], color=INKB, lw=1.1,
                            solid_capstyle="round")

        # column heads, row frames, arrows and row names
        HEAD = ("Mendelian consistency in trios", "Copy number in trios", "Copy number in field monitoring")
        for i, (x, h) in enumerate(zip(CX, HEAD)):
            t(x - 0.3, 107.8, "abc"[i], size=S.LETTER, fontweight="bold")
            t(x + 4.0, 107.8, h, size=T1, fontweight="bold")
            frame(x, YD[0], CW, YD[1] - YD[0]); frame(x, YE[0], CW, YE[1] - YE[0])
            arr(x + CW / 2, YD[0], x + CW / 2, YC[1]); arr(x + CW / 2, YC[0], x + CW / 2, YE[1])
        frame(CX[0], YC[0], CX[2] + CW - CX[0], YC[1] - YC[0])                 # all five callers, run on every part
        for (y0, y1), s in ((YD, "Data"), (YC, "Callers"), (YE, "Evaluation")):
            t(2.8, (y0 + y1) / 2, s, size=T1, ha="center", rotation=90, fontweight="bold")
        TY = YD[1] - 3.6                                                         # title line inside a data box

        # b, data: the four H. zea trios (square male, circle female) with their Tryp77 copies by ddPCR
        x0 = CX[1]
        t(x0 + 2.2, TY, f"Four {S.ZEA} trios", size=T2, fontweight="bold")
        for (tr, cns), (x, y) in zip(TRIOS.items(), ((x0 + 8.5, 92.5), (x0 + 33.0, 92.5),
                                                     (x0 + 8.5, 79.5), (x0 + 33.0, 79.5))):
            trio(x, y, cns, SEX[tr])
        copies(x0 + 2.2, 69.9, 1); t(x0 + 4.4, 70.3, f"1 {S.TRYP} copy (ddPCR)")

        # a, data: each trio moth (here a child) is called with the nine moths of the other three trios as its panel;
        # its own parents (dashed) are never in the panel
        x0 = CX[0]
        t(x0 + 2.2, TY, "Moth called against a panel", size=T2, fontweight="bold")
        frame(x0 + 14.0, 76.0, 36.0, 19.5, fc="#ececec", ec="none")
        for i, tr in enumerate(SEX):
            x, y, dx, dy = x0 + 5.0 + i * 11.6, 91.5, 6.4, 9.5
            q = i == 0
            ax.plot([x, x + dx], [y, y], color=FAINT if q else INKB, lw=S.THIN + 0.1)
            ax.plot([x + dx / 2] * 2, [y, y - dy], color=FAINT if q else INKB, lw=S.THIN + 0.1)
            for (xx, yy), sx, who in zip(((x, y), (x + dx, y), (x + dx / 2, y - dy)), SEX[tr], "ppc"):
                if q and who == "c":
                    person(xx, yy, sx, r=1.3, fc=INKB)
                elif q:
                    person(xx, yy, sx, r=1.3, ec="#8c8c8c", ls=(0, (1.2, 0.8)))
                else:
                    person(xx, yy, sx, r=1.3)
        t(x0 + 8.2, 72.0, "Moth", ha="center")
        t(x0 + 32.0, 72.0, "Panel: other trios", ha="center")

        # c, data: the Louisiana moths by year, and the panel growing as each year is added
        x0 = CX[2]
        t(x0 + 2.2, TY, "Louisiana field moths", size=T2, fontweight="bold")
        YRS = (("2002", 6, 0), ("2012", 6, 1), ("2017", 6, 4))                 # cartoon: year, moths, carriers
        SP = 2.3                                                                # dot spacing
        YTOP = tuple(94.5 - i * (SP + 4.4) for i in range(3))                   # top row of each year
        PCX = (x0 + 18.0, x0 + 31.0, x0 + 44.0)                                 # panel centres
        for i, (yr, n, c) in enumerate(YRS):
            t(x0 + 2.2, YTOP[i] - SP * (math.ceil(n / 3) - 1) / 2, yr)
        for p, cx in enumerate(PCX):
            bot = YTOP[p] - SP * (math.ceil(YRS[p][1] / 3) - 1) - 1.4
            frame(cx - 4.4, bot, 8.8, YTOP[0] + 1.5 - bot, fc="#ececec", ec="none")
            for i, (yr, n, c) in enumerate(YRS[:p + 1]):
                for k in range(n):
                    r_, q = divmod(k, 3)
                    ax.scatter(cx - SP + q * SP, YTOP[i] - r_ * SP, s=7, fc=INKB if k < c else "white", ec=INKB,
                               lw=0.5, zorder=3)
            if p:
                arr(PCX[p - 1] + 4.6, YTOP[0] - 1.2, cx - 4.6, YTOP[0] - 1.2)
        ax.scatter([x0 + 3.0], [70.3], s=7, fc=INKB, ec=INKB, lw=0.5)
        t(x0 + 4.4, 70.3, f"{S.TRYP} carrier")
        t(x0 + CW - 2.2, 70.3, "Growing panel", ha="right")

        # callers: one box, as every part runs all of them; two halves drawn alike (title, pictogram, callers)
        X = lambda x0, w: (lambda x: x0 + (x - 61) * w / 45)                    # sketches drawn on 61..106, rescaled
        HX = (CX[0], (CX[0] + CX[2] + CW) / 2)                                  # left edges of the two halves
        PW, YP = 48.0, 48.5                                                     # pictogram width and baseline
        ax.plot([HX[1] - 1.5] * 2, [YC[0] + 2.0, YC[1] - 2.0], color=FAINT, lw=S.THIN)
        for hx, title, sub, ms in ((HX[0], "Breakpoint callers", "split and discordant reads", ("Delly2", "smoove")),
                                   (HX[1], "Read-depth callers", "depth against a panel of other moths",
                                    ("GATK-gCNV", "CNVcaller", "CNVkit"))):
            t(hx + 2.2, YC[1] - 3.4, title, size=T2, fontweight="bold")
            t(hx + 2.2, YC[1] - 7.2, sub)
            for k, m in enumerate(ms):
                yy = (YC[0] + YC[1]) / 2 + (len(ms) - 1) / 2 * 4.6 - k * 4.6
                ax.add_patch(mp.Circle((hx + PW + 12.0, yy), 1.1, fc=S.COLORS[m], ec=INKB, lw=0.4))
                t(hx + PW + 13.9, yy, m, size=8.0)
        tt = X(HX[0] + 4.0, PW)                                                 # split and discordant reads over a deletion
        y = YP
        ax.plot([tt(62), tt(76)], [y, y], color=INKB, lw=1.3); ax.plot([tt(88), tt(104)], [y, y], color=INKB, lw=1.3)
        ax.plot([tt(76), tt(88)], [y, y], color=INKB, lw=0.6, ls=(0, (1, 1.2)))
        ax.add_patch(mp.Rectangle((tt(76), y - 1.2), tt(88) - tt(76), 2.4, fc="#ececec", ec="none"))
        for (a_, b_), yy in ((((70, 76), (88, 93)), y + 2.0), (((65, 70.5), (92, 97.5)), y + 4.0)):
            for xa, xb in (a_, b_):
                ax.plot([tt(xa), tt(xb)], [yy, yy], color=MID, lw=1.6, solid_capstyle="butt")
            ax.plot([tt(a_[1]), tt(b_[0])], [yy, yy], color=MID, lw=0.5, ls=(0, (2, 1.5)))
        BIAS = lambda x: 0.55 * np.sin((x - 61) * 0.9) + 0.35 * np.sin((x - 61) * 2.3)  # depth bias shared by all moths
        tt = X(HX[1] + 4.0, PW)                                                 # the moth's depth over the panel's range
        xs = np.linspace(61, 106, 90); f = (xs - 61) / 45
        ax.fill_between([tt(v) for v in xs], y - 1.3 + 0.6 * BIAS(xs), y + 1.3 + 0.6 * BIAS(xs), color="#d9d9d9", lw=0)
        ax.step([tt(v) for v in xs], y + 0.6 * BIAS(xs) + np.where((f > 0.45) & (f < 0.65), 3.4, 0)
                + np.random.default_rng(7).normal(0, 0.12, xs.size), where="mid", color=INKB, lw=0.8)

        ET = YE[1] - 3.6                                                        # title line inside an evaluation box
        # b, evaluation: called copy number at Tryp77 against ddPCR
        x0 = CX[1]
        t(x0 + 2.2, ET, f"Copy number at {S.TRYP}", size=T2, fontweight="bold")
        callers(x0 + CW - 2.4, ET, CN_CALLERS)
        px, py, pw, ph = x0 + 14.0, 7.5, 30.0, 21.0
        ax.plot([px, px, px + pw], [py + ph, py, py], color=INKB, lw=S.LW)
        ax.plot([px, px + pw], [py, py + ph], color=INKB, lw=S.THIN, ls=(0, (1.5, 1.5)))
        for a, b in ((0.12, 0.16), (0.16, 0.10), (0.20, 0.22), (0.34, 0.28), (0.46, 0.52), (0.58, 0.50),
                     (0.70, 0.74), (0.80, 0.70), (0.90, 0.86)):                 # cartoon points, not data
            ax.scatter(px + a * pw, py + b * ph, s=10, fc="#a6a6a6", ec=INKB, lw=0.4, zorder=3)
        t(px + pw / 2, py - 3.4, "ddPCR", ha="center")
        t(px - 3.0, py + ph / 2, "Caller", ha="center", rotation=90)

        # a, evaluation: is each child call found in a parent?
        x0 = CX[0]
        t(x0 + 2.2, ET, "Child calls in parents", size=T2, fontweight="bold")
        callers(x0 + CW - 2.4, ET, ORDER)
        tr = {"P1": 26.0, "P2": 19.5, "C": 13.0}
        for k, sex, lab in (("P1", "M", "Parent"), ("P2", "F", "Parent"), ("C", "F", "Child")):
            ax.plot([x0 + 15.0, x0 + 47.0], [tr[k], tr[k]], color="#a6a6a6", lw=S.THIN)
            person(x0 + 12.3, tr[k], sex, r=1.1)
            t(x0 + 2.2, tr[k], lab)
        for k, a, b in (("P1", 18.5, 26.5), ("P2", 29.0, 33.5), ("C", 19.5, 26.5), ("C", 37.5, 42.5)):
            ax.add_patch(mp.Rectangle((x0 + a, tr[k] - 0.9), b - a, 1.8, fc=MID, ec="none"))
        ax.plot([x0 + 23.0] * 2, [tr["C"] + 0.9, tr["P1"] - 0.9], color=INKB, lw=S.THIN, ls=(0, (1, 1)))
        tick(x0 + 23.0, 6.5, True); tick(x0 + 40.0, 6.5, False)

        # c, evaluation: carriers called in each panel (true carriers outlined) and chr 7 sites rising over the years
        x0 = CX[2]
        t(x0 + 2.2, ET, "Carriers and rising sites", size=T2, fontweight="bold")
        callers(x0 + CW - 2.4, ET, MON)
        bx, by, bh = x0 + 8.0, 7.5, 21.0
        ax.plot([bx - 1.5, bx - 1.5, bx + 19.5], [by + bh, by, by], color=INKB, lw=S.LW)
        for k, (yr, f) in enumerate((("2002", 0.0), ("2012", 0.2), ("2017", 0.85))):   # true carriers, cartoon
            xb = bx + 1.0 + k * 6.5
            ax.add_patch(mp.Rectangle((xb, by), 4.0, max(f * bh, 0.3), fc="none", ec=INKB, lw=S.THIN + 0.1,
                                      ls=(0, (2, 1))))
            t(xb + 2.0, by - 3.0, yr, ha="center")
        t(bx - 4.2, by + bh / 2, "Carriers", ha="center", rotation=90)
        fx, fw = x0 + 33.0, 17.0
        ax.plot([fx, fx, fx + fw], [by + bh, by, by], color=INKB, lw=S.LW)
        for k, (a, b) in enumerate(((6.0, 6.6), (10.0, 9.4), (14.0, 14.4))):  # sites whose frequency does not change
            ax.plot([fx + 1, fx + fw / 2, fx + fw - 1], [by + a, by + (a + b) / 2 + 0.6 * (-1) ** k, by + b],
                    color="#a6a6a6", lw=S.THIN + 0.1)
        ax.plot([fx + 1, fx + fw / 2, fx + fw - 1], [by + 2.0, by + 3.5, by + 19.0], color=INKB, lw=1.1)
        t(fx + fw - 2.6, by + 19.5, S.TRYP, ha="right", color=INKB)
        t(fx + 1, by - 3.0, "2002", ha="left"); t(fx + fw, by - 3.0, "2017", ha="right")
        t(fx - 2.6, by + bh / 2, "Frequency", ha="center", rotation=90)
        S.save(fig, "study_design")


# ── 2. scorecard ──────────────────────────────────────────────────────────────
def caller_scorecard():
    MIR = {(r["species"], r["caller"]): float(r["mir"]) for r in tsv("mir_headline")
           if r["svtype"] == "DUP" and r["calls"] == "filtered"}
    DDR = {r["caller"]: float(r["pearson_r"]) for r in tsv("tryp77_ddpcr_correlation") if r["caller"] != "smoove"}
    TRUTH = {r["caller"]: r for r in tsv("truthset_human") if (r["truth"], r["calls"], r["svtype"]) == ("HGSVC2", "raw", "DEL")}
    # false carriers: false duplications plus false deletions among the 35 moths of the full time series (Additional file 1, Fig. S26)
    FALSE = {m: int(PC["2002+2012+2017", MON_ROW[m]]["false_dup"]) + int(PC["2002+2012+2017", MON_ROW[m]]["false_del"])
             for m in MON}
    # set-up time on chromosome 7 (hours, median per panel; Table S20): GATK-gCNV model training, CNVcaller discovery
    # plus genotyping, Delly2 merge-genotype-classify of a panel.
    RT = {(r["caller"], r["step"]): r["wall_med_h"] for r in runtime.summarise(runtime.load())
          if r["species"] == "Monitoring"}
    SETUP = {"GATK-gCNV": RT["GATK-gCNV", "model training, COHORT (chromosome 7, max CN 10)"],
             "CNVcaller": RT["CNVcaller", "CNV-region discovery (time-series and first constant-size panels)"]
                          + RT["CNVcaller", "genotyping (time-series and first constant-size panels)"],
             "Delly2": RT["Delly2", "delly merge, genotyping (delly cnv), classify; all panels"]}
    ADDS = {"GATK-gCNV": 1, "CNVcaller": 0, "Delly2": 1}       # a new moth is called without re-running the panel
    ZEA, HUM, MONG, COST = r"$H.\ zea$ trios", "Human trios", "Monitoring", "Running it"
    # (group, header, values, better, kind, scale max, format); kind: bar | dots1 | dots8 | cross | logbar | check
    COLS = [
        (ZEA, "Dup. MIR", {m: MIR["H. zea", m] for m in ORDER}, "low", "bar", 50, "{:.1f}%"),
        (ZEA, "ddPCR $r$", DDR, "high", "bar", 1, "{:.2f}"),       # smoove reports no copy number
        (HUM, "Dup. MIR", {m: MIR["Human", m] for m in ORDER}, "low", "bar", 50, "{:.1f}%"),
        (HUM, "Del.\nprecision", {m: float(TRUTH[m]["precision"]) for m in ORDER}, "high", "bar", 100, "{:.0f}%"),
        (HUM, "Del.\nrecall", {m: float(TRUTH[m]["recall"]) for m in ORDER}, "high", "bar", 100, "{:.0f}%"),
        (MONG, "rare", RARE, "high", "dots1", None, None),
        (MONG, "common", COMMON, "high", "dots8", None, None),
        (MONG, "False\ncarriers", FALSE, "low", "cross", None, None),
        (MONG, "Sites\nflagged", FLAGGED, "low", "bar", max(FLAGGED.values()), "{:.0f}"),
        (COST, "Set-up", SETUP, "low", "logbar", None, None),
        (COST, "Adds a moth\nwithout re-run", ADDS, "high", "check", None, None),
    ]
    W, H = 170, 62
    fig, ax = canvas(W, H)
    X0, CW, RH = 27.5, 11.4, 9.0
    XB = X0 + len(COLS) * CW + 0.8                                                     # "best for" column
    ry = {m: 42.5 - i * RH - (2.0 if i >= 2 else 0) for i, m in enumerate(ORDER)}     # row centres
    for i, m in enumerate(ORDER):
        if i % 2 == 0:
            box(ax, 0.5, ry[m] - RH / 2 + 0.2, XB + 13 - 0.5, RH - 0.4, fc="#F6F6F6", r=1.2)


    def fmt_h(h):
        return f"{h * 60:.0f} min" if h < 1 else f"{h:g} h"


    for j, (grp, head, vals, better, kind, top, fmt) in enumerate(COLS):
        cx, xl, Lw = X0 + (j + 0.5) * CW, X0 + j * CW + 1.4, CW - 2.8
        best = (min if better == "low" else max)(vals.values())
        for m in ORDER:
            v, cy = vals.get(m), ry[m]
            if v is None:
                continue
            bold = v == best
            if kind in ("bar", "logbar"):
                f = (math.log10(v) + 2) / 3 if kind == "logbar" else v / top               # logbar: 0.6 min .. 10 h
                ax.add_patch(mp.Rectangle((xl, cy + 0.3), Lw, 2.0, fc="#E6E6E6", ec="none"))
                ax.add_patch(mp.Rectangle((xl, cy + 0.3), max(Lw * min(f, 1), 0.35), 2.0, fc=COLORS[m], ec="none"))
                s_ = fmt_h(v) if kind == "logbar" else fmt.format(v)
                txt(ax, cx, cy - 2.0, s_, size=5.2, ha="center", color=INK if bold else MUTED,
                    fontweight="bold" if bold else "normal")
            elif kind == "dots1":
                dot(ax, cx, cy, v, r=1.5, color=COLORS[m])
            elif kind == "dots8":
                for k in range(8):
                    dot(ax, cx - 3.3 + (k % 4) * 2.2, cy + (1.1 if k < 4 else -1.1), k < v, r=0.8, color=COLORS[m])
            elif kind == "cross":
                if v == 0:
                    txt(ax, cx, cy, "0", size=6.2, ha="center", fontweight="bold")
                for k in range(v):
                    txt(ax, cx + (k - (v - 1) / 2) * 2.6, cy, "✗", size=8, ha="center", color=COLORS[m], **SYM)
            elif kind == "check":
                txt(ax, cx, cy, "✓" if v else "✗", size=8.5, ha="center", color=INK if v else MUTED, **SYM)
        txt(ax, cx, 48.6, head, size=5.6, ha="center", va="bottom")
        if kind not in ("dots1", "dots8", "check"):
            txt(ax, cx, 47.3, "↓" if better == "low" else "↑", size=5.0, ha="center", color=MUTED, **SYM)
    # one span header over the two carrier columns
    j0 = next(j for j, c_ in enumerate(COLS) if c_[4] == "dots1")
    txt(ax, X0 + (j0 + 1) * CW, 51.6, r"$Tryp77$ carriers called", size=5.6, ha="center", va="bottom")
    ax.plot([X0 + j0 * CW + 1.4, X0 + (j0 + 2) * CW - 1.4], [51.2, 51.2], color=MUTED, lw=0.4)

    groups = []
    for j, c_ in enumerate(COLS):
        if not groups or groups[-1][0] != c_[0]:
            groups.append([c_[0], j, j])
        groups[-1][2] = j
    for g, a, b in groups:
        xa, xb = X0 + a * CW + 0.5, X0 + (b + 1) * CW - 0.5
        ax.plot([xa, xb], [56.4, 56.4], color=INK, lw=0.5)
        txt(ax, xa + 0.5, 58.6, g, size=6.6, fontweight="bold")
        if a:
            ax.plot([xa - 0.5] * 2, [0.6, 55.0], color="#DDDDDD", lw=0.5)

    # caller names with a paradigm pictogram (split read / depth step)
    for m in ORDER:
        y = ry[m]
        ax.add_patch(mp.Circle((2.2, y), 1.0, fc=COLORS[m], ec="none"))
        txt(ax, 4.0, y, m, size=6.6)
        if m in ("Delly2", "smoove"):
            ax.plot([20.0, 22.1], [y, y], color=DARK, lw=1.4, solid_capstyle="butt")
            ax.plot([24.1, 26.2], [y, y], color=DARK, lw=1.4, solid_capstyle="butt")
            ax.plot([22.1, 24.1], [y, y], color=DARK, lw=0.4, ls=(0, (1, 0.8)))
        else:
            ax.plot([20.0, 21.8, 21.8, 24.4, 24.4, 26.2], [y - 0.8, y - 0.8, y + 0.9, y + 0.9, y - 0.8, y - 0.8], color=DARK, lw=0.7)

    # best monitoring use
    txt(ax, XB + 6, 48.6, "Best for", size=5.6, ha="center", va="bottom")
    for m, use in (("GATK-gCNV", "rare CNV"), ("CNVcaller", "common CNV"), ("Delly2", "unknown CNV")):
        txt(ax, XB + 6, ry[m], use, size=5.8, ha="center")
    save(fig, "caller_scorecard", max_mm=None)


if __name__ == "__main__":
    study_design()
    caller_scorecard()
