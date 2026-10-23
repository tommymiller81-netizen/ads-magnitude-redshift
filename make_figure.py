"""
Reproduce the Hubble vs anti-de Sitter comparison figure for the DESI DR1
bright-time sample.

Input : desi_dr1_main_bright_sample.csv.gz, written by desi_bright_selection.py
Output: desi_bright_hubble_vs_ads.png (300 dpi) and .pdf

Usage:  python make_figure.py [path/to/desi_dr1_main_bright_sample.csv.gz]
Requires: numpy, pandas, matplotlib.
"""
import sys
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

mpl.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.titlelocation": "left", "lines.linewidth": 1.2, "patch.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42})

f3  = lambda x: format(x, "#.3g").rstrip(".")      # 3 significant figures
neg = lambda s: s.replace("-", "\u2212")           # typographic minus

def main(path="desi_dr1_main_bright_sample.csv.gz"):
    d = pd.read_csv(path)
    z = d.z.to_numpy(float); m = 22.5 - 2.5*np.log10(d.flux_r.to_numpy(float))
    lz = np.log10(z); isQ = (d.spectype == "QSO").to_numpy()
    M = {"M_H": m - 5*lz - 43.0, "M_AdS": m - 2.5*np.log10(1 - (1 + z)**-2) - 55.3}

    edges = np.arange(-2.6, 0.75, 0.1); cen = 0.5*(edges[1:] + edges[:-1]); idx = np.digitize(lz, edges) - 1
    bandcols = ["#fffbe6", "#fff3bf", "#ffe88a"]; edgec = "#c9a227"
    cG, cQ, cF, dimc = "#0072B2", "#D55E00", "#000000", "#333333"     # Okabe-Ito blue/vermillion
    xarr = {"99.7": -2.54, "95": -2.45, "68": -2.36}; L = 0.9
    ap = dict(arrowstyle="-|>", lw=0.6, color=dimc, shrinkA=0, shrinkB=0, mutation_scale=5)
    med = {k: float(np.median(v)) for k, v in M.items()}; offset = med["M_AdS"] - med["M_H"]

    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.6), constrained_layout=True)
    for ax, (key, title) in zip(axs, [("M_H", "Hubble (q\u2080 = 1, C = \u221243.0)"),
                                      ("M_AdS", "Anti-de Sitter (eq. 21, C = \u221255.3)")]):
        v = M[key]; p = np.percentile(v, [0.15, 2.5, 16, 50, 84, 97.5, 99.85])
        for (lo, hi), c in zip([(0, 6), (1, 5), (2, 4)], bandcols):
            ax.axhspan(p[lo], p[hi], color=c, lw=0, zorder=0)
        for y in (p[0], p[1], p[2], p[4], p[5], p[6]):
            ax.axhline(y, color=edgec, lw=0.4, alpha=0.6, zorder=1)
        hb = ax.hexbin(lz, v, gridsize=(110, 70), extent=(-2.6, 0.75, p[3] - 9, p[3] + 9), bins="log",
                       cmap="Greys", mincnt=3, linewidths=0, zorder=2)
        sl, ic = np.polyfit(lz, v, 1); xx = np.array([-2.6, 0.75])
        ax.plot(xx, ic + sl*xx, color=cF, lw=1.1, ls=(0, (5, 2)), zorder=4)
        def bmed(sel, nmin):
            vv, ii = v[sel], idx[sel]
            return np.array([np.median(vv[ii == b]) if (ii == b).sum() >= nmin else np.nan for b in range(len(cen))])
        ax.plot(cen, bmed(isQ, 200), "s-", color=cQ, ms=3.2, mfc="none", mew=0.8, lw=0.9, zorder=5)
        ax.plot(cen, bmed(~isQ, 1000), "o-", color=cG, ms=2.6, lw=1.1, zorder=6)
        for lab, (lo, hi) in [("68", (2, 4)), ("95", (1, 5)), ("99.7", (0, 6))]:
            x, ytop, ybot = xarr[lab], p[lo], p[hi]
            ax.annotate("", xy=(x, ytop), xytext=(x, ytop - L), arrowprops=ap, zorder=8)
            ax.annotate("", xy=(x, ybot), xytext=(x, ybot + L), arrowprops=ap, zorder=8)
            ax.annotate(f"{f3(ybot - ytop)} mag", xy=(x, ytop - L), xytext=(0, 1.5), textcoords="offset points",
                        rotation=90, ha="center", va="bottom", fontsize=5.5, color=dimc, zorder=8)
        yc = med["M_H"] if key == "M_H" else med["M_H"] + offset
        ax.set_ylim(yc + 9, yc - 9); ax.set_xlim(-2.6, 0.75); ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.text(0.97, 0.03, f"least-squares slope (all objects) = {neg(f3(sl))} mag per dex",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6, zorder=7)
        ax.set_title(title); ax.set_xlabel(r"$\log_{10} z$"); ax.set_ylabel(r"absolute magnitude $M_r$")
    cb = fig.colorbar(hb, ax=axs, shrink=0.75, pad=0.01); cb.set_label("objects per cell")
    h = [Patch(facecolor=bandcols[2], edgecolor=edgec, lw=0.4, label="central 68%"),
         Line2D([], [], color=cG, marker="o", ms=2.6, lw=1.1, label="galaxy median per 0.1 dex"),
         Patch(facecolor=bandcols[1], edgecolor=edgec, lw=0.4, label="central 95%"),
         Line2D([], [], color=cQ, marker="s", ms=3.2, mfc="none", mew=0.8, lw=0.9, label="quasar median per 0.1 dex"),
         Patch(facecolor=bandcols[0], edgecolor=edgec, lw=0.4, label="central 99.7%"),
         Line2D([], [], color=cF, lw=1.1, ls=(0, (5, 2)), label="least-squares fit, all objects")]
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.47, 0.0), ncol=3, fontsize=6, columnspacing=1.6)
    for ax, letter in zip(axs, "ab"):
        ax.text(-0.18, 1.02, letter, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom", ha="left")
    nG, nQ = int((~isQ).sum()), int(isQ.sum())
    fig.text(0.01, -0.115,
        f"DESI DR1 main, bright time (DESI Collaboration 2025, arXiv:2503.14745): N = {len(d):,} ({nG:,} galaxies, {nQ:,} quasars) with ZCAT_PRIMARY and ZWARN = 0;\n"
        "Legacy Surveys r, no K- or extinction correction.\n"
        "All objects retained; the thin vertical streaks near log z \u2248 0.01 and 0.20\u20130.23 are redshift outliers. Median bins need \u2265 1000 galaxies or \u2265 200 quasars.\n"
        f"Both y-axes span 18 mag; panel b is offset from panel a by {neg(f3(offset))} mag. Arrowheads mark the edges of each central band; labels give its full width.",
        fontsize=6, va="top")
    for ext in ("png", "pdf"):
        fig.savefig(f"desi_bright_hubble_vs_ads.{ext}")
    print("wrote desi_bright_hubble_vs_ads.png/.pdf")

if __name__ == "__main__":
    main(*sys.argv[1:2])
