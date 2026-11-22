"""Print style of the main-text figures: font, sizes, line widths, caller palette, panel letters and saving.

Figures are drawn at their printed size (Genome Biology: 170 mm full width, 85 mm half width, at most
225 mm high), so the sizes below are the printed sizes. Text is Times New Roman bold, as in the RawMod
paper figures: boxed axes with inward ticks on all four sides, a light horizontal grid, black-edged
markers and bars, and framed legends. The smallest text anywhere is SMALL (7.5 pt); panel letters are
lower-case bold at FS + 3.

A figure function draws inside `with figstyle.paper():`, so the style applies to that figure only and
the other figures of the same script keep the style they were drawn with. Fonts: Times New Roman TTFs
in CNVB_FONTS or ~/.fonts are used when present; otherwise a metric-compatible serif."""
import glob
import os
from contextlib import contextmanager
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.lines
import matplotlib.pyplot as plt
import matplotlib.text
from matplotlib.colors import to_rgb

from .paths import FIGURES

MM = 1 / 25.4
FULL, HALF, MAX_H = 170.0, 85.0, 225.0          # printed widths and maximum height (mm)

for _d in (os.environ.get("CNVB_FONTS", ""), os.path.expanduser("~/.fonts")):
    for _f in glob.glob(os.path.join(_d, "*.[tT][tT][fF]")) if _d else []:
        fm.fontManager.addfont(_f)
FONT = next((n for n in ("Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif")
             if any(f.name == n for f in fm.fontManager.ttflist)), "serif")

FS = 9.0                # axis labels and panel titles (pt, printed)
TICK = 8.5              # tick labels
SMALL = 7.5             # legends, in-panel annotations: the smallest text used
LETTER = FS + 3         # panel letters
LW = 0.8                # axes, ticks and reference lines (pt)
THIN = 0.5              # marker edges, grid

STYLE = {
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "font.family": FONT, "font.weight": "bold",
    "axes.labelweight": "bold", "axes.titleweight": "bold",
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS,
    "xtick.labelsize": TICK, "ytick.labelsize": TICK, "legend.fontsize": SMALL, "legend.title_fontsize": SMALL,
    "mathtext.fontset": "custom", "mathtext.rm": f"{FONT}:bold", "mathtext.it": f"{FONT}:bold:italic",
    "mathtext.bf": f"{FONT}:bold", "mathtext.default": "bf",
    "mathtext.sf": f"{FONT}:bold:italic",       # italic() below: unlike \\it, keeps digits italic (Tryp77)
    "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
    "axes.edgecolor": "black", "axes.linewidth": LW,
    "axes.spines.top": True, "axes.spines.right": True,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.major.width": LW, "ytick.major.width": LW, "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.minor.width": THIN, "ytick.minor.width": THIN, "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
    "xtick.major.pad": 3, "ytick.major.pad": 3, "axes.labelpad": 2.5, "axes.titlepad": 3,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": "#b0b0b0", "grid.linestyle": "-",
    "grid.linewidth": 0.4, "axes.axisbelow": True,
    "lines.linewidth": 1.2, "lines.markeredgewidth": THIN, "patch.linewidth": THIN,
    "errorbar.capsize": 1.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "legend.frameon": True, "legend.framealpha": 1.0, "legend.edgecolor": "black", "legend.fancybox": False,
    "legend.handlelength": 1.6, "legend.handletextpad": 0.4, "legend.columnspacing": 1.2,
    "legend.borderpad": 0.35, "legend.labelspacing": 0.3, "legend.borderaxespad": 0.3,
}

# One colour per caller, the RawMod / RawHash2 hues in each caller's earlier hue family (orange, red,
# purple, blue, green). Breakpoint callers are drawn dashed where lines are used.
COLORS = {"Delly2": "#db8043", "smoove": "#c0392b", "GATK-gCNV": "#aa4499", "CNVcaller": "#4b71bb",
          "CNVkit": "#1baf7a"}
ORDER = ["Delly2", "smoove", "GATK-gCNV", "CNVcaller", "CNVkit"]      # breakpoint | read depth
BREAKPOINT = ["Delly2", "smoove"]
DASH = (0, (3.2, 1.6))
GREY, LIGHT = "#7f7f7f", "#d9d9d9"
EDGE = dict(edgecolor="black", linewidth=THIN)


def italic(s):
    """Text set in bold italic as mathtext, digits included (gene and species names): italic("Tryp77")."""
    return r"$\mathsf{" + s.replace(" ", r"\ ") + "}$"


TRYP = italic("Tryp77")
ZEA = italic("H. zea")


def linestyle(caller):
    return DASH if caller in BREAKPOINT else "-"


def tint(color, f):
    """Blend a colour toward white; f = 1 is the colour itself."""
    r, g, b = to_rgb(color)
    return (1 - f + f * r, 1 - f + f * g, 1 - f + f * b)


@contextmanager
def paper():
    """Draw (and save) a figure in the print style."""
    with plt.rc_context(STYLE):
        yield


def label(ax, s, x=-0.12, y=1.0, va="bottom", ha="left", **kw):
    """Panel letter: lower-case bold, no dot or parentheses, placed in axes fraction (default: on the
    axes' top edge, to the left of it)."""
    return ax.text(x, y, s, transform=ax.transAxes, fontsize=LETTER, fontweight="bold", va=va, ha=ha, **kw)


def save(fig, name, width_mm=FULL, max_h=MAX_H):
    """Write results/figures/{name}.pdf and .png. The canvas is width_mm wide and saved without
    bbox cropping, so the PDF prints at exactly that width; the height must stay within max_h."""
    w, h = fig.get_size_inches() / MM
    assert abs(w - width_mm) < 0.5, f"{name}: drawn {w:.1f} mm wide, expected {width_mm} mm"
    assert h <= max_h + 0.5, f"{name}: {h:.0f} mm high (max {max_h} mm)"
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())      # nothing may fall outside the canvas
    W, H = fig.get_size_inches()
    out = [bb.x0 < -0.01, bb.y0 < -0.01, bb.x1 > W + 0.01, bb.y1 > H + 0.01]
    if any(out):
        print(f"warning: {name}: artists outside the canvas (left, bottom, right, top) = {out}; "
              f"bbox {bb.x0 / MM:.1f}, {bb.y0 / MM:.1f}, {bb.x1 / MM:.1f}, {bb.y1 / MM:.1f} mm")
    sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    widths = [a.get_linewidth() for a in fig.findobj(matplotlib.lines.Line2D) if a.get_visible() and a.get_linewidth()]
    if min(sizes) < SMALL:
        print(f"warning: {name}: text at {min(sizes):g} pt (smallest allowed {SMALL} pt)")
    if widths and min(widths) < 0.25:
        print(f"warning: {name}: line of {min(widths):g} pt (< 0.25 pt)")
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{name}.pdf")
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    plt.close(fig)
    print("wrote", FIGURES / f"{name}.pdf", f"({w:.0f} x {h:.0f} mm; smallest text {min(sizes):g} pt)")
