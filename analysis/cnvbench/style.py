"""Figure style of the supplementary figures: sizes, panel letters, saving. The caller palette comes from
figstyle, the print style of the main-text figures (which supplementary figures can adopt later)."""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

from .figstyle import COLORS, ORDER                                  # one caller palette for every figure
from .paths import FIGURES

MM = 1 / 25.4
INK, MUTED, GREY = "#222222", "#6B6B6B", "#4D4D4D"

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 7, "axes.titlesize": 7.5,
    "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "grid.linewidth": 0.4,
    "grid.color": "#DDDDDD", "axes.edgecolor": "#555555", "text.color": INK,
    "axes.labelcolor": INK, "xtick.color": "#444444", "ytick.color": "#444444",
})
# Times New Roman when its TTFs are installed (the paper's figures); matplotlib's serif otherwise.
for f in glob.glob(os.path.expanduser("~/.fonts/[Tt]imes*.[Tt][Tt][Ff]")):
    fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
                     "mathtext.fontset": "custom", "mathtext.rm": "serif", "mathtext.it": "serif:italic",
                     "mathtext.bf": "serif:bold"})


def tint(color, f):
    """Blend a colour toward white; f = 1 is the colour itself."""
    r, g, b = to_rgb(color)
    return (1 - f + f * r, 1 - f + f * g, 1 - f + f * b)


def label(ax, s, x=-0.2, y=1.04):
    """Bold panel letter."""
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def fit_width(fig, max_mm=169.5, pad_in=0.1):
    """Narrow the canvas until the tight bounding box fits max_mm (Genome Biology: 170 mm page width)."""
    for _ in range(8):
        fig.canvas.draw()
        w = (fig.get_tightbbox(fig.canvas.get_renderer()).width + 2 * pad_in) * 25.4
        if w <= max_mm:
            return
        fig.set_size_inches(fig.get_size_inches()[0] - (w - max_mm + 0.3) / 25.4, fig.get_size_inches()[1])


def save(fig, name, max_mm=169.5):
    """Write results/figures/{name}.pdf and .png at no more than max_mm wide."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    if max_mm:
        fit_width(fig, max_mm)
    for ext in ("pdf", "png"):
        fig.savefig(FIGURES / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", FIGURES / f"{name}.pdf")
