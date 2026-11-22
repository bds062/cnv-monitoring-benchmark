"""CNVcaller discovery frequency -f 0.1 (the paper's setting) against -f 0: does the same moth's Tryp77 call
change as the panel around it changes? Every monitoring CNVcaller run was repeated with only -f changed
('f0' stage of workflow/monitoring/run_all.sh), into the monitoring folder of a second results tree,
CNVB_RESULTS_F0 (default: results_f0 next to CNVB_RESULTS). Additional file 1, Section S4.4 and Fig. S15; tables cnvcaller_f0_calls (every
call at both settings) and cnvcaller_f0_summary (carriers called, false calls and call changes per context).
Truth: ddPCR for the trio moths, read depth for the published moths (carrier = copy number >= 2.5)."""
import os
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

import cnvbench.monitoring as M
from cnvbench import figstyle as F
from cnvbench.locus import ddpcr
from cnvbench.paths import RESULTS

ROOT = {"0.1": M.MON, "0": Path(os.environ.get("CNVB_RESULTS_F0", RESULTS.parent / "results_f0")) / "monitoring"}
SETTINGS = ("0.1", "0")
DD = ddpcr()
DEP = M.read_depth()
CARRIER = {s: cn >= M.DUP for s, (_, cn) in DEP.items()}
CARRIER.update({s: DD[s] >= M.DUP for _, _, s in M.TRIO_MOTHS})
TRIO = {s: (t, r) for t, r, s in M.TRIO_MOTHS}
K = (0, 2, 4, 6, 8)
SERIES = range(1, 11)
NEST = [f"nest_s{r}_k{k}" for k in K for r in SERIES]
GROWING = ("2002", "2002+2012", "2002+2012+2017")
SWEEP = sorted((p for p in M.PANELS if p.startswith("sweep_")), key=lambda p: tuple(map(int, re.findall(r"\d+", p))))
OTHER = [p for p in M.PANELS if not p.startswith(("nest_", "sweep_"))]
state = lambda v: "dup" if v >= M.DUP else "del" if v <= M.DEL else "none"


def calls(root):
    """{(kind, panel, moth): copy number}: every panel member, and every trio moth queried against a panel."""
    out = {}
    for p in NEST + SWEEP + OTHER:
        for s, v in M.cnvcaller_cp(root / "panels" / p / "cnvcaller" / f"{p}.vcf.gz").items():
            out["panel", p, s] = v
    for p in NEST + list(GROWING):
        for s in TRIO:
            out["query", p, s] = M.cnvcaller_cp(root / "query" / p / "cnvcaller" / f"{s}.vcf.gz")[s]
    return out


CN = {f: calls(ROOT[f]) for f in SETTINGS}
assert CN["0.1"].keys() == CN["0"].keys()
M.write_tsv("cnvcaller_f0_calls", ("kind", "panel", "moth", "carrier", "cn_f0.1", "cn_f0"),
            [(*key, int(CARRIER[key[2]]), CN["0.1"][key], CN["0"][key]) for key in CN["0"]])

# contexts of Fig. S15 c-d (rows top to bottom); call changes along the dilution series (table only)
CONTEXTS = [("Queried, 0/11 panels", lambda k, p: k == "query" and p.endswith("_k0")),
            ("Queried, 2–8/11 panels", lambda k, p: k == "query" and p.startswith("nest") and not p.endswith("_k0")),
            ("Queried, growing panels", lambda k, p: k == "query" and p in GROWING),
            ("Members, dilution panels", lambda k, p: k == "panel" and p.startswith("nest")),
            ("Members, time-series panels", lambda k, p: k == "panel" and p in OTHER)]
CONTEXTS += [(f"Members, sweep {n}/13 ", lambda k, p, n=n: k == "panel" and p.startswith(f"sweep_k{n}_")) for n in (1, 2, 3, 4)]


def tally(f, sel):
    n = Counter(M.outcome(CARRIER[s], v) for (k, p, s), v in CN[f].items() if sel(k, p))
    return n["tp"], n["tp"] + n["fn"], n["fp"] + n["del"], n["tn"] + n["fp"] + n["del"]


def changes(f, kind, carrier):
    """(changed, compared): same moth in consecutive nested panels of a series (k -> k + 2)."""
    ch = tot = 0
    for r in SERIES:
        for k in K[:-1]:
            a, b = f"nest_s{r}_k{k}", f"nest_s{r}_k{k + 2}"
            for s in {s for (kk, p, s) in CN[f] if kk == kind and p == a and CARRIER[s] == carrier}:
                if (kind, b, s) in CN[f]:
                    tot += 1
                    ch += state(CN[f][kind, a, s]) != state(CN[f][kind, b, s])
    return ch, tot


CHANGE_ROWS = [("Queried carriers", "query", True), ("Queried non-carriers", "query", False),
               ("Member carriers", "panel", True), ("Member non-carriers", "panel", False)]
rows = []
for name, sel in CONTEXTS:
    for f in SETTINGS:
        rows.append((name, f, *tally(f, sel), "", ""))
for name, kind, car in CHANGE_ROWS:
    for f in SETTINGS:
        rows.append(("Call changes, " + name, f, "", "", "", "", *changes(f, kind, car)))
M.write_tsv("cnvcaller_f0_summary", ("context", "f", "carriers_called", "carriers", "false_calls", "non_carriers",
            "calls_changed", "calls_compared"), rows)
for r in rows:
    print("\t".join(map(str, r)))
flips = {f: sum(len({state(v) for (k, p, s2), v in CN[f].items() if s2 == s}) > 1 for s in CARRIER) for f in SETTINGS}
print("moths whose call is not the same in every panel: " + ", ".join(f"-f {f}: {n}/{len(CARRIER)}" for f, n in flips.items()))

# figure (print style of the paper, laid out in mm)
OK_DARK, OK_LIGHT, BAD_DARK, BAD_LIGHT = "#006B4F", "#C6EADC", "#D55E00", "#EE9A6A"     # as Fig. 4 and Fig. S13
CODE = {"tp": 0, "tn": 1, "fn": 2, "fp": 3, "del": 3}
CMAP = ListedColormap([OK_DARK, OK_LIGHT, BAD_DARK, BAD_LIGHT])
SET_COL = {"0.1": F.tint(F.COLORS["CNVcaller"], 0.4), "0": F.COLORS["CNVcaller"]}
SET_MARK = {"0.1": "o", "0": "s"}
SET_LAB = {"0.1": "CNVcaller -f 0.1 (as in the paper)", "0": "CNVcaller -f 0"}
ROLE = {"resistant": "res. parent", "susceptible": "sus. parent", "child": "child"}
ORDER_M = sorted(TRIO, key=lambda s: -DD[s])
GAP = 1.6
XPOS = [i * len(SERIES) + j + i * GAP for i in range(len(K)) for j in range(len(SERIES))]

W, L_, R_ = F.FULL, 40.0, 1.0                    # mm; L_ holds the row names
CELL, TITLE, TICKS, XLAB = 3.0, 5.0, 4.6, 5.0     # heatmap row; title, tick-label and axis-label strips
ROW = 4.4                                        # one context row of c-d
LEG = 11.0                                       # legend strip
HEAT = len(ORDER_M) * CELL
H = 0.5 + 2 * (TITLE + HEAT) + 1.5 + TICKS + XLAB + 3.0 + TITLE + len(CONTEXTS) * ROW + TICKS + XLAB + 2.0 + LEG + 0.5


def axes(fig, x, y, w, h):
    """Axes at (x, y, w, h) in mm from the bottom left."""
    return fig.add_axes([x / W, y / H, w / W, h / H])


def title(fig, x, y, text, letter, lx=None):
    fig.text(x / W, (y + 1.0) / H, text, ha="left", va="bottom", fontsize=F.FS, fontweight="bold")
    fig.text((lx if lx is not None else 1.0) / W, (y + 1.0) / H, letter, ha="left", va="bottom",
             fontsize=F.LETTER, fontweight="bold")


def heatmap(fig, y0, f, letter, xlabel):
    """The 12 trio moths (rows, by ddPCR copy number) queried against the 50 nested panels (columns)."""
    ax = axes(fig, L_, y0, W - L_ - R_, HEAT)
    grid = np.array([[CODE[M.outcome(CARRIER[s], CN[f]["query", p, s])] for p in NEST] for s in ORDER_M])
    for c, x in enumerate(XPOS):
        ax.imshow(grid[:, [c]], cmap=CMAP, vmin=-0.5, vmax=3.5, aspect="auto", interpolation="nearest",
                  extent=(x, x + 0.9, len(ORDER_M), 0))
    for k in range(len(K)):                      # frame each block of ten series
        x = XPOS[k * len(SERIES)]
        ax.add_patch(plt.Rectangle((x - 0.05, 0), len(SERIES) - 0.0, len(ORDER_M), fill=False, lw=F.LW,
                                   edgecolor="black", zorder=4))
    n_car = sum(CARRIER[s] for s in ORDER_M)
    for k in range(len(K)):                      # carriers above the line, non-carriers below
        x = XPOS[k * len(SERIES)]
        ax.plot([x - 0.05, x + len(SERIES) - 0.05], [n_car, n_car], color="black", lw=F.LW, zorder=4)
    ax.set_xlim(-0.4, XPOS[-1] + 1.3); ax.set_ylim(len(ORDER_M), 0)
    ax.set_yticks(np.arange(len(ORDER_M)) + 0.5,
                  [f"{TRIO[s][0].replace('Trio', 'Trio ')} {ROLE[TRIO[s][1]]} ({DD[s]:.1f})" for s in ORDER_M])
    ax.set_xticks([XPOS[j * len(SERIES)] + len(SERIES) / 2 for j in range(len(K))], [f"{k}/11" for k in K] if xlabel else [])
    ax.tick_params(length=0, top=False, right=False)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    if xlabel:
        ax.set_xlabel("Carriers in the 11-moth carrier-dilution panel (series 1–10 in each block)")
    title(fig, L_, y0 + HEAT, f"{SET_LAB[f]}: trio moths queried against carrier-dilution panels", letter)


def dots(fig, x0, y0, w, num, den, ttl, letter, lx, names):
    """Carriers called (num, den = 0, 1) or non-carriers falsely called (2, 3) per context, both settings."""
    ax = axes(fig, x0, y0, w, len(CONTEXTS) * ROW)
    y = np.arange(len(CONTEXTS))[::-1]
    t = {f: [tally(f, sel) for _, sel in CONTEXTS] for f in SETTINGS}
    for yy, a, b in zip(y, t["0.1"], t["0"]):
        if a[den]:
            ax.plot([100 * a[num] / a[den], 100 * b[num] / b[den]], [yy, yy], color=F.GREY, lw=F.LW, zorder=2)
    for f in SETTINGS:
        v = [100 * a[num] / a[den] if a[den] else np.nan for a in t[f]]
        ax.scatter(v, y, s=26, marker=SET_MARK[f], color=SET_COL[f], zorder=3, **F.EDGE)
    for yy, a in zip(y, t["0"]):
        ax.text(1.03, yy, f"{a[num]}/{a[den]}", fontsize=F.SMALL, va="center", transform=ax.get_yaxis_transform())
    ax.set_yticks(y, [n for n, _ in CONTEXTS] if names else [])
    ax.tick_params(axis="y", left=False, right=False)
    ax.grid(axis="x"); ax.grid(axis="y", visible=False)
    ax.set_xlim((-6, 106) if num == 0 else (-2, 30))
    ax.set_ylim(-0.6, len(CONTEXTS) - 0.4)
    ax.set_xlabel("Calls (%)")
    title(fig, x0, y0 + len(CONTEXTS) * ROW, ttl, letter, lx)


with F.paper():
    fig = plt.figure(figsize=(W * F.MM, H * F.MM))
    y = H - 0.5 - TITLE - HEAT
    heatmap(fig, y, "0.1", "a", False)
    y -= 1.5 + TITLE + HEAT
    heatmap(fig, y, "0", "b", True)
    y -= TICKS + XLAB + 3.0 + TITLE + len(CONTEXTS) * ROW
    dots(fig, L_, y, 52, 0, 1, "Carriers called", "c", None, True)
    dots(fig, L_ + 52 + 17, y, 44, 2, 3, "Non-carriers falsely called", "d", L_ + 52 + 11, False)
    sw = lambda c, m="s", e="black": plt.Line2D([], [], ls="", marker=m, ms=6, mfc=c, mec=e, mew=F.THIN)
    hs = [(sw(OK_DARK), "Carrier, called"), (sw(BAD_DARK), "Carrier, missed"),
          (sw(OK_LIGHT), "Non-carrier, not called"), (sw(BAD_LIGHT), "Non-carrier, false dup. or del."),
          (sw(SET_COL["0.1"], "o"), SET_LAB["0.1"]), (sw(SET_COL["0"], "s"), SET_LAB["0"])]
    fig.legend([h for h, _ in hs], [t for _, t in hs], loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.5 / H),
               handlelength=1.2, columnspacing=1.4)
    F.save(fig, "cnvcaller_f0_flipping")
