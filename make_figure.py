"""
Reproduce the figures of the paper.

1. Main figure: absolute magnitude against log10 z. Rows: (a) SDSS, (b) DESI DR1 bright time.
   Columns: (1) anti-de Sitter, (2) Hubble (q0 = 1).
   Output: absolute_magnitude_sdss_desi.pdf/.png
2. Appendix figure: apparent magnitude m_r against log10 z, (a) SDSS and (b) DESI.
   Output: apparent_magnitude_sdss_desi.pdf/.png
Inputs: desi_dr1_main_bright_sample.csv.gz (written by desi_bright_selection.py) and
        sdss_legacy_sample.csv.gz (written by sdss_legacy_selection.py).

Panel letters name the sample (a = SDSS, b = DESI) in both figures; digits name the relation.
Gray shading is the fraction of each sample per cell, so samples of different size share one
scale. Every panel has the same physical scale (1.405 dex/in, 6.80 mag/in), so slopes can be
compared by eye across all panels. Use the PDFs for JHEP; the PNGs are previews. The figure
captions are in the paper (LaTeX), not in the images.

Usage:  python make_figure.py [desi.csv.gz] [sdss.csv.gz]
Requires: numpy, pandas, matplotlib.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm
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
sgn = lambda x: neg(f"{x:+.1f}")                                 # one decimal, explicit sign, typographic minus


BANDC = ["#fffbe6", "#fff3bf", "#ffe88a"]; EDGEC = "#dfc77d"        # opaque (no transparency layers)
CG, CQ, CF, DIMC = "#0072B2", "#D55E00", "#000000", "#7a5c0f"     # Okabe-Ito blue/vermillion; deep gold arrows
XARR = {"99.7": -2.40, "95": -2.14, "68": -1.88}; L = 1.1           # band-width arrows (log z positions, length in mag)
SPAN = 9.0                                                          # each panel shows median +/- 9 mag

FIT_DASH_PT = (4.0, 2.0)   # user's choice (option C)
# dash pattern of the least-squares line, in units of its line width (1.1 pt): FIT_DASH_PT dash / gap in points.
# The legend sample is 2.0 x 8 pt = 16 pt long, which the pattern fills with whole dashes (3 x 4 + 2 x 2).
FIT_DASH = tuple(v/1.1 for v in FIT_DASH_PT)

def _panel(ax, P, ylabel, invert, ylim, xarr, xarr_lo, label_at, ap, edges, cen, first, xlabel=True, weights=None):
    """Draw one panel (bands, hexbin, fit, medians, arrows, limits, titles); returns (hexbin, (slope, intercept)).
    weights: per-object weights for the hexbin (e.g. 1/N, the fraction of the sample per cell); None = counts."""
    hbs, fits = [], []
    lz, v, isQ = P["lz"], P["v"], P["isQ"]; idx = np.digitize(lz, edges) - 1
    p = np.percentile(v, [0.15, 2.5, 16, 50, 84, 97.5, 99.85])
    for (lo, hi), c in zip([(0, 6), (1, 5), (2, 4)], BANDC):
        ax.axhspan(p[lo], p[hi], color=c, lw=0, zorder=0)
    for y in (p[0], p[1], p[2], p[4], p[5], p[6]):
        ax.axhline(y, color=EDGEC, lw=0.4, zorder=1)
    if P.get("ylim") or ylim:
        y0, y1 = P.get("ylim") or ylim                              # fixed range for this panel / all panels
        P["ycentre"] = 0.5*(y0 + y1)
    else:
        # median +/- SPAN, shifted only if needed so the outer arrows (and the width label on the
        # top arrow) stay inside the panel
        c = p[3]; room_lab = 0.9; pad = 0.2
        if invert:   # brighter (smaller) values at the top, fainter at the bottom
            c = min(c, p[0] - L - room_lab + SPAN)   # top arrow and its label fit:  c - SPAN <= p0 - L - room
            c = max(c, p[6] + L + pad - SPAN)        # bottom arrow fits:            c + SPAN >= p6 + L + pad
        y0, y1 = c - SPAN, c + SPAN
        P["ycentre"] = c
    ny = int(round(70*(y1 - y0)/(2*SPAN)))                             # same cell size as an 18-mag panel
    if weights is None:
        hb = ax.hexbin(lz, v, gridsize=(110, ny), extent=(-2.6, 0.75, y0, y1), bins="log",
                       cmap="Greys", mincnt=3, linewidths=0, zorder=2)
    else:   # summed weights per cell; mincnt still counts objects (cells with >= 3 objects are shown)
        hb = ax.hexbin(lz, v, C=weights, reduce_C_function=np.sum, gridsize=(110, ny),
                       extent=(-2.6, 0.75, y0, y1), norm=LogNorm(), cmap="Greys", mincnt=3,
                       linewidths=0, zorder=2)
    hbs.append(hb)
    sl, ic = np.polyfit(lz, v, 1); xx = np.array([-2.6, 0.75]); fits.append((sl, ic))
    ax.plot(xx, ic + sl*xx, color=CF, lw=1.1, ls=(0, FIT_DASH), zorder=4)
    def bmed(sel, nmin):
        vv, ii = v[sel], idx[sel]
        return np.array([np.median(vv[ii == b]) if (ii == b).sum() >= nmin else np.nan for b in range(len(cen))])
    medG, medQ = bmed(~isQ, NMIN_G), bmed(isQ, NMIN_Q)
    ax.plot(cen, medQ, "-", color=CQ, lw=0.9, zorder=9)
    ax.plot(cen, medQ, ls="none", marker="s", mfc=CQ, mec=CQ, mew=0.3, ms=3.0, zorder=9)
    ax.plot(cen, medG, "o-", color=CG, ms=2.6, lw=1.1, zorder=10)  # medians above arrows (zorder 8)
    for lab, (lo, hi) in [("68", (2, 4)), ("95", (1, 5)), ("99.7", (0, 6))]:
        # e_up / e_dn: band edges that appear at the top / bottom of the panel
        e_up, e_dn = (p[lo], p[hi]) if invert else (p[hi], p[lo])
        s_ = -1 if invert else 1                                   # +1 if larger values plot higher
        xu, xd = xarr[lab], xarr_lo[lab]
        ax.annotate("", xy=(xu, e_up), xytext=(xu, e_up + s_*L), arrowprops=ap, zorder=8)   # points down
        ax.annotate("", xy=(xd, e_dn), xytext=(xd, e_dn - s_*L), arrowprops=ap, zorder=8)   # points up
        # width label above the upper arrow, or below the lower arrow
        if label_at == "top":
            ax.annotate(f"{p[hi] - p[lo]:.1f}", xy=(xu, e_up + s_*L), xytext=(0, 1.5), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8, color=DIMC, zorder=8)
        else:
            ax.annotate(f"{p[hi] - p[lo]:.1f}", xy=(xd, e_dn - s_*L), xytext=(0, -1.5), textcoords="offset points",
                        ha="center", va="top", fontsize=8, color=DIMC, zorder=8)
    ax.set_ylim(*((y1, y0) if invert else (y0, y1))); ax.set_xlim(-2.6, 0.75); ax.yaxis.set_major_locator(MultipleLocator(2))
    if P.get("title"): ax.set_title(P["title"])
    if xlabel: ax.set_xlabel(r"$\log_{10} z$")
    else: ax.tick_params(labelbottom=False)
    if first: ax.set_ylabel(ylabel)
    return hbs[0], fits[0]


def draw_grid(rows, out, ylabels, letters, invert=True, xarr=None, ymag=18.0, target_w=2.384):
    """Panels in a grid: rows[r][c] is a panel dict (lz, v, isQ, ylim, title or None). Each panel is
    target_w wide and ymag/6.8 in tall, so all figures share 1.405 dex/in and 6.80 mag/in.
    ylabels[r] labels row r (on its first panel); letters[r][c] is the panel letter.
    Gray shading: fraction of the panel's sample per cell (weights 1/N), one color scale for all panels.
    Legend and color bar share one strip below the panels, no wider than the panels."""
    ap = dict(arrowstyle="-|>,head_length=0.5,head_width=0.25", lw=1.4, color=DIMC, shrinkA=0, shrinkB=0, mutation_scale=10)
    xarr = xarr or XARR
    edges = np.arange(-2.6, 0.75, 0.1); cen = 0.5*(edges[1:] + edges[:-1])
    nr, nc = len(rows), len(rows[0]); target = (target_w, ymag/6.8)
    fig, axs = plt.subplots(nr, nc, figsize=(5.85, 3.0*nr), constrained_layout=True, squeeze=False)
    hbs = []
    for r, row in enumerate(rows):
        for c, P in enumerate(row):
            ax = axs[r, c]; w = np.full(len(P["v"]), 1.0/len(P["v"]))
            hb, fit = _panel(ax, dict(P), ylabels[r], invert, None, xarr, xarr, "top", ap, edges, cen, c == 0,
                             xlabel=(r == nr - 1), weights=w)
            hbs.append(hb)
            top = ax.get_ylim()[1]                                  # value at the upper edge (inverted or not)
            yl = top + XL_TOP[1] if invert else top - XL_TOP[1]
            ax.text(0.5, yl, f"slope = {sgn(fit[0])}", transform=ax.get_yaxis_transform(),
                    ha="center", va="bottom", fontsize=8, zorder=7)
            ax.text(-0.18, 1.02, letters[r][c], transform=ax.transAxes, fontweight="bold", fontsize=10,
                    va="bottom", ha="left")
    h = [Patch(facecolor=BANDC[2], edgecolor=EDGEC, lw=0.4, label="central 68%"),
         Line2D([], [], color=CG, marker="o", ms=2.6, lw=1.1, label="galaxy median"),
         Patch(facecolor=BANDC[1], edgecolor=EDGEC, lw=0.4, label="central 95%"),
         Line2D([], [], color=CQ, marker="s", ms=3.0, mfc=CQ, mec=CQ, mew=0.3, lw=0.9, label="quasar median"),
         Patch(facecolor=BANDC[0], edgecolor=EDGEC, lw=0.4, label="central 99.7%"),
         Line2D([], [], color=CF, lw=1.1, ls=(0, FIT_DASH), label="least-squares fit")]
    # ncol=3 fills column by column: row 1 = bands, row 2 = galaxy median, quasar median, LSQ fit
    leg = fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.0, 0.0), ncol=3, fontsize=8, columnspacing=0.9,
                     handletextpad=0.5, borderaxespad=0.0)
    lo = min(h_.get_clim()[0] for h_ in hbs); hi = max(h_.get_clim()[1] for h_ in hbs)
    for h_ in hbs: h_.set_norm(LogNorm(lo, hi))
    # panels at the common printed size: adjust the figure size until the axes match the target
    for _ in range(6):
        fig.canvas.draw(); b_ = axs[0, 0].get_position(); W, H = fig.get_size_inches()
        fig.set_size_inches(W + nc*(target[0] - b_.width*W), H + nr*(target[1] - b_.height*H))
    # legend and color bar side by side; the color-bar label ends at the right edge of the panels
    fig.canvas.draw(); W, H = fig.get_size_inches()
    x_ax1 = axs[0, nc - 1].get_position().x1*W
    gap, cbw, cbh = 0.25, 1.15, 0.086
    bb = leg.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    cy = bb.y0 + 0.62*bb.height
    cax = fig.add_axes([0.5, cy/H, cbw/W, cbh/H])
    cb = fig.colorbar(hbs[0], cax=cax, orientation="horizontal"); cb.set_label("fraction of sample per cell", fontsize=8)
    cb.set_ticks([1e-6, 1e-5, 1e-4, 1e-3]); cb.ax.tick_params(labelsize=7)
    fig.canvas.draw()
    lab = cb.ax.xaxis.label.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    over = max(0.0, 0.5*(lab.width - cbw))
    xb = x_ax1 - over - cbw
    cax.set_position([xb/W, cy/H, cbw/W, cbh/H])
    leg.set_bbox_to_anchor(((xb - over - gap - bb.width)/W, 0.0), transform=fig.transFigure)
    # solid frame outlines (matplotlib's default "no outline" is written as a transparent state)
    fig.patch.set_edgecolor(fig.patch.get_facecolor())
    for a in fig.axes: a.patch.set_edgecolor(a.patch.get_facecolor())
    fig.canvas.draw(); W, H = fig.get_size_inches(); b_ = axs[0, 0].get_position()
    print(f"{out}: axes {b_.width*W:.3f} x {b_.height*H:.3f} in; {3.35/(b_.width*W):.3f} dex/in, "
          f"{ymag/(b_.height*H):.3f} mag/in; color scale {lo:.2e} to {hi:.2e}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}")
    plt.close(fig); print(f"wrote {out}.png/.pdf")


def load_desi(path):
    d = pd.read_csv(path)
    z = d.z.to_numpy(float); m = 22.5 - 2.5*np.log10(d.flux_r.to_numpy(float))
    return z, m, (d.spectype == "QSO").to_numpy()

def load_sdss(path):
    s = pd.read_csv(path)
    return s.z_use.to_numpy(float), s.modelMag_r.to_numpy(float), (s["sample"] == "QSO").to_numpy()

def main(desi="desi_dr1_main_bright_sample.csv.gz", sdss="sdss_legacy_sample.csv.gz"):
    z, m, isQ = load_desi(desi); lz = np.log10(z)
    zs, ms, isQs = load_sdss(sdss); lzs = np.log10(zs)
    # y-ranges, 18 mag each: Hubble centred near the DESI median (-20.6), anti-de Sitter offset by exactly
    # 14 mag (a multiple of the 2-mag tick spacing) so the ticks of both columns line up
    yl_H = (-29.5, -11.5); yl_A = (yl_H[0] - 14.0, yl_H[1] - 14.0)
    tA, tH = r"Anti-de Sitter ($C_{\rm AdS} = -55.3$)", r"Hubble ($C_{\rm H} = -43.0$)"
    ads = lambda z_, m_: m_ - 2.5*np.log10(1 - (1 + z_)**-2) - 55.3
    hub = lambda z_, m_: m_ - 5*np.log10(z_) - 43.0
    rows = [[dict(lz=lzs, v=ads(zs, ms), isQ=isQs, ylim=yl_A, title=tA),
             dict(lz=lzs, v=hub(zs, ms), isQ=isQs, ylim=yl_H, title=tH)],
            [dict(lz=lz, v=ads(z, m), isQ=isQ, ylim=yl_A, title=None),     # column titles on the top row only
             dict(lz=lz, v=hub(z, m), isQ=isQ, ylim=yl_H, title=None)]]
    draw_grid(rows, "absolute_magnitude_sdss_desi",
              [r"$\mathbf{SDSS}$ absolute magnitude", r"$\mathbf{DESI}$ absolute magnitude"],
              [["a1", "a2"], ["b1", "b2"]])
    # appendix: apparent magnitude, one row, (a) SDSS and (b) DESI; not inverted (fainter upward);
    # 15-mag range at the same mag/in, so slopes compare directly with the main figure
    draw_grid([[dict(lz=lzs, v=ms, isQ=isQs, ylim=(11.0, 26.0), title="SDSS"),
                dict(lz=lz, v=m, isQ=isQ, ylim=(11.0, 26.0), title="DESI")]],
              "apparent_magnitude_sdss_desi", ["Apparent magnitude"], [["a", "b"]],
              invert=False, xarr=XARR_APP, ymag=15.0)

XL_TOP = ("top", 2.8)    # slope labels in all figures: horizontal, centred, bottom 2.8 mag (0.41 in) below the top edge
# appendix arrows: opposing arrows share one vertical line; slightly wider spacing than the main
# figure because the DESI band edges bunch up at the faint limit
XARR_APP = {"99.7": -2.45, "95": -2.08, "68": -1.71}

if __name__ == "__main__":
    main(*sys.argv[1:3])
