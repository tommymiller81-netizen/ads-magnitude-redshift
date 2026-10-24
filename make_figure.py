"""
Reproduce the Hubble vs anti-de Sitter comparison figure for the DESI DR1
bright-time sample.

Input : desi_dr1_main_bright_sample.csv.gz, written by desi_bright_selection.py
Output: desi_bright_hubble_vs_ads.pdf (vector, fonts embedded; use this for JHEP) and a .png preview
The figure caption is in the paper (LaTeX), not in the image.

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
from matplotlib.path import Path

# Hollow square drawn as a filled ring, so the PDF needs no transparency (JHEP asks for none)
_o = [(-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1)]
_i = [(0.62*x, 0.62*y) for x, y in _o[::-1]]
RING = Path(_o + _i, [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY] + [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY])

mpl.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.titlelocation": "left", "lines.linewidth": 1.2, "patch.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42})

f3  = lambda x: format(x, "#.3g").rstrip(".")      # 3 significant figures
NMIN_G, NMIN_Q = 1000, 200                          # min. galaxies / quasars per bin for running medians
neg = lambda s: s.replace("-", "\u2212")           # typographic minus

def main(path="desi_dr1_main_bright_sample.csv.gz"):
    d = pd.read_csv(path)
    z = d.z.to_numpy(float); m = 22.5 - 2.5*np.log10(d.flux_r.to_numpy(float))
    lz = np.log10(z); isQ = (d.spectype == "QSO").to_numpy()
    M = {"M_H": m - 5*lz - 43.0, "M_AdS": m - 2.5*np.log10(1 - (1 + z)**-2) - 55.3}

    edges = np.arange(-2.6, 0.75, 0.1); cen = 0.5*(edges[1:] + edges[:-1]); idx = np.digitize(lz, edges) - 1
    bandcols = ["#fffbe6", "#fff3bf", "#ffe88a"]; edgec = "#dfc77d"   # opaque (no transparency layers)
    cG, cQ, cF, dimc = "#0072B2", "#D55E00", "#000000", "#7a5c0f"     # Okabe-Ito blue/vermillion; deep gold arrows (band hue, >=5:1 contrast on white/bands)
    xarr = {"99.7": -2.40, "95": -2.14, "68": -1.88}; L = 1.1
    ap = dict(arrowstyle="-|>,head_length=0.5,head_width=0.25", lw=1.4, color=dimc, shrinkA=0, shrinkB=0, mutation_scale=10)
    med = {k: float(np.median(v)) for k, v in M.items()}; offset = med["M_AdS"] - med["M_H"]

    # drawn at the printed size: JHEP \textwidth = 0.72 x A4 width = 151 mm = 5.95 in, so include at width=\textwidth (scale ~1)
    fig, axs = plt.subplots(1, 2, figsize=(5.85, 3.3), constrained_layout=True)
    hbs = []
    fits = {}
    for ax, (key, title) in zip(axs, [("M_AdS", r"Anti-de Sitter ($C_{\rm AdS} = -55.3$)"),
                                      ("M_H", r"Hubble ($C_{\rm H} = -43.0$)")]):
        v = M[key]; p = np.percentile(v, [0.15, 2.5, 16, 50, 84, 97.5, 99.85])
        for (lo, hi), c in zip([(0, 6), (1, 5), (2, 4)], bandcols):
            ax.axhspan(p[lo], p[hi], color=c, lw=0, zorder=0)
        for y in (p[0], p[1], p[2], p[4], p[5], p[6]):
            ax.axhline(y, color=edgec, lw=0.4, zorder=1)
        hb = ax.hexbin(lz, v, gridsize=(110, 70), extent=(-2.6, 0.75, p[3] - 9, p[3] + 9), bins="log",
                       cmap="Greys", mincnt=3, linewidths=0, zorder=2); hbs.append(hb)
        sl, ic = np.polyfit(lz, v, 1); xx = np.array([-2.6, 0.75])
        ax.plot(xx, ic + sl*xx, color=cF, lw=1.1, ls=(0, (5, 2)), zorder=4)
        def bmed(sel, nmin):
            vv, ii = v[sel], idx[sel]
            return np.array([np.median(vv[ii == b]) if (ii == b).sum() >= nmin else np.nan for b in range(len(cen))])
        medG, medQ = bmed(~isQ, NMIN_G), bmed(isQ, NMIN_Q)
        ax.plot(cen, medQ, "-", color=cQ, lw=0.9, zorder=9)
        ax.plot(cen, medQ, ls="none", marker="s", mfc=cQ, mec=cQ, mew=0.3, ms=3.0, zorder=9)
        ax.plot(cen, medG, "o-", color=cG, ms=2.6, lw=1.1, zorder=10)  # medians above arrows (zorder 8)
        for lab, (lo, hi) in [("68", (2, 4)), ("95", (1, 5)), ("99.7", (0, 6))]:
            x, ytop, ybot = xarr[lab], p[lo], p[hi]
            ax.annotate("", xy=(x, ytop), xytext=(x, ytop - L), arrowprops=ap, zorder=8)
            ax.annotate("", xy=(x, ybot), xytext=(x, ybot + L), arrowprops=ap, zorder=8)
            ax.annotate(f"{ybot - ytop:.1f}", xy=(x, ytop - L), xytext=(0, 1.5), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8, color=dimc, zorder=8)
        yc = med["M_H"] if key == "M_H" else med["M_H"] + offset
        ax.set_ylim(yc + 9, yc - 9); ax.set_xlim(-2.6, 0.75); ax.yaxis.set_major_locator(MultipleLocator(2))
        fits[key] = (sl, ic)
        ax.set_title(title); ax.set_xlabel(r"$\log_{10} z$")
        if ax is axs[0]: ax.set_ylabel("Absolute magnitude")
    h = [Patch(facecolor=bandcols[2], edgecolor=edgec, lw=0.4, label="central 68%"),
         Line2D([], [], color=cG, marker="o", ms=2.6, lw=1.1, label="galaxy median"),
         Patch(facecolor=bandcols[1], edgecolor=edgec, lw=0.4, label="central 95%"),
         Line2D([], [], color=cQ, marker="s", ms=3.0, mfc=cQ, mec=cQ, mew=0.3, lw=0.9, label="quasar median"),
         Patch(facecolor=bandcols[0], edgecolor=edgec, lw=0.4, label="central 99.7%"),
         Line2D([], [], color=cF, lw=1.1, ls=(0, (5, 2)), label="least-squares fit")]
        # ncol=3 fills column by column: row 1 = bands, row 2 = galaxy median, quasar median, LSQ fit
    fig.legend(handles=h, loc="upper center",
               bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=8, columnspacing=1.6)
    # horizontal color bar below the legend, outside the panel layout so the panels keep full width
    cax = fig.add_axes([0.33, -0.205, 0.30, 0.026])
    # one color scale for both panels (identical cell size: same extent and gridsize)
    lo = min(h_.get_clim()[0] for h_ in hbs); hi = max(h_.get_clim()[1] for h_ in hbs)
    for h_ in hbs: h_.set_clim(lo, hi)
    cb = fig.colorbar(hbs[0], cax=cax, orientation="horizontal"); cb.set_label("objects per cell", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    for ax, letter in zip(axs, "ab"):
        ax.text(-0.18, 1.02, letter, transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom", ha="left")
    # slope labels: parallel to each fit line, just above it, near the sparse left end.
    # The angle is computed after the layout is final so it matches the drawn line.
    fig.canvas.draw()
    for ax, key in zip(axs, ("M_AdS", "M_H")):
        sl, ic = fits[key]; xl = {"M_H": -1.40, "M_AdS": -1.98}[key]; yl = ic + sl*xl
        p1, p2 = ax.transData.transform([(xl - 0.1, ic + sl*(xl - 0.1)), (xl + 0.1, ic + sl*(xl + 0.1))])
        ang = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])
        ax.annotate(f"slope\u2002{neg(f'{sl:.1f}')}", xy=(xl, yl),
                    xytext=(-4*np.sin(ang), 4*np.cos(ang)), textcoords="offset points",
                    rotation=np.degrees(ang), rotation_mode="anchor", ha="center", va="bottom",
                    fontsize=8, zorder=7)
    # solid frame outlines (matplotlib's default "no outline" is written as a transparent state)
    fig.patch.set_edgecolor(fig.patch.get_facecolor())
    for a in fig.axes: a.patch.set_edgecolor(a.patch.get_facecolor())   # panels and color bar
    for ext in ("png", "pdf"):
        fig.savefig(f"desi_bright_hubble_vs_ads.{ext}")
    print("wrote desi_bright_hubble_vs_ads.png/.pdf")

if __name__ == "__main__":
    main(*sys.argv[1:2])
