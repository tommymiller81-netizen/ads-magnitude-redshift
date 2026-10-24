"""
Reproduce the figures of the paper.

1. Main figure: absolute magnitude against log10 z for the DESI DR1 bright-time sample,
   (a) anti-de Sitter and (b) Hubble (q0 = 1).
   Input : desi_dr1_main_bright_sample.csv.gz (written by desi_bright_selection.py)
   Output: desi_bright_hubble_vs_ads.pdf/.png
2. Appendix figure: apparent magnitude m_r against log10 z, (a) SDSS and (b) DESI.
   Input : the DESI file above and sdss_legacy_sample.csv.gz (written by sdss_legacy_selection.py)
   Output: apparent_magnitude_sdss_desi.pdf/.png  (made only if the SDSS file exists)

Every panel spans 18 mag vertically at the same physical size, so slopes can be compared
by eye across all panels of both figures. Use the PDFs for JHEP; the PNGs are previews.
The figure captions are in the paper (LaTeX), not in the images.

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

def draw(panels, ylabel, out, invert=True, ylim=None, xarr=None, xarr_lo=None, figsize=(5.85, 3.3), cbar_y=-0.205,
         label_at="top"):
    """panels: list of dicts with keys lz, v, isQ, title, xl (log z of the slope label).
    invert=True puts brighter (smaller) magnitudes at the top, the convention for absolute magnitude;
    invert=False is the classic Hubble-diagram orientation for apparent magnitude (fainter upward)."""
    ap = dict(arrowstyle="-|>,head_length=0.5,head_width=0.25", lw=1.4, color=DIMC, shrinkA=0, shrinkB=0, mutation_scale=10)
    xarr = xarr or XARR; xarr_lo = xarr_lo or xarr   # x of the arrows at the band edges drawn at the top / bottom of the panel
    edges = np.arange(-2.6, 0.75, 0.1); cen = 0.5*(edges[1:] + edges[:-1])
    # drawn at the printed size: JHEP \textwidth = 0.72 x A4 width = 151 mm = 5.95 in, so include at width=\textwidth (scale ~1)
    fig, axs = plt.subplots(1, len(panels), figsize=figsize, constrained_layout=True)
    hbs, fits = [], []
    for ax, P in zip(axs, panels):
        lz, v, isQ = P["lz"], P["v"], P["isQ"]; idx = np.digitize(lz, edges) - 1
        p = np.percentile(v, [0.15, 2.5, 16, 50, 84, 97.5, 99.85])
        for (lo, hi), c in zip([(0, 6), (1, 5), (2, 4)], BANDC):
            ax.axhspan(p[lo], p[hi], color=c, lw=0, zorder=0)
        for y in (p[0], p[1], p[2], p[4], p[5], p[6]):
            ax.axhline(y, color=EDGEC, lw=0.4, zorder=1)
        y0, y1 = ylim if ylim else (p[3] - SPAN, p[3] + SPAN)           # common range, or median +/- SPAN
        ny = int(round(70*(y1 - y0)/(2*SPAN)))                             # same cell size as an 18-mag panel
        hb = ax.hexbin(lz, v, gridsize=(110, ny), extent=(-2.6, 0.75, y0, y1), bins="log",
                       cmap="Greys", mincnt=3, linewidths=0, zorder=2); hbs.append(hb)
        sl, ic = np.polyfit(lz, v, 1); xx = np.array([-2.6, 0.75]); fits.append((sl, ic))
        ax.plot(xx, ic + sl*xx, color=CF, lw=1.1, ls=(0, (5, 2)), zorder=4)
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
        ax.set_title(P["title"]); ax.set_xlabel(r"$\log_{10} z$")
        if ax is axs[0]: ax.set_ylabel(ylabel)
    h = [Patch(facecolor=BANDC[2], edgecolor=EDGEC, lw=0.4, label="central 68%"),
         Line2D([], [], color=CG, marker="o", ms=2.6, lw=1.1, label="galaxy median"),
         Patch(facecolor=BANDC[1], edgecolor=EDGEC, lw=0.4, label="central 95%"),
         Line2D([], [], color=CQ, marker="s", ms=3.0, mfc=CQ, mec=CQ, mew=0.3, lw=0.9, label="quasar median"),
         Patch(facecolor=BANDC[0], edgecolor=EDGEC, lw=0.4, label="central 99.7%"),
         Line2D([], [], color=CF, lw=1.1, ls=(0, (5, 2)), label="least-squares fit")]
    # ncol=3 fills column by column: row 1 = bands, row 2 = galaxy median, quasar median, LSQ fit
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=8, columnspacing=1.6)
    # horizontal color bar below the legend, outside the panel layout so the panels keep full width
    cax = fig.add_axes([0.33, cbar_y, 0.30, 0.026*3.3/figsize[1]])
    # one color scale for all panels of a figure (identical cell size: same extent and gridsize)
    lo = min(h_.get_clim()[0] for h_ in hbs); hi = max(h_.get_clim()[1] for h_ in hbs)
    for h_ in hbs: h_.set_clim(lo, hi)
    cb = fig.colorbar(hbs[0], cax=cax, orientation="horizontal"); cb.set_label("objects per cell", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    for ax, letter in zip(axs, "ab"):
        ax.text(-0.18, 1.02, letter, transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom", ha="left")
    # slope labels: parallel to each fit line, just above it. The angle is computed after the
    # layout is final so it matches the drawn line.
    fig.canvas.draw()
    for ax, P, (sl, ic) in zip(axs, panels, fits):
        if isinstance(P["xl"], tuple):                           # ("centre", y): horizontal label, centred, bottom at y
            ax.text(0.5, P["xl"][1], f"slope = {sgn(sl)}", transform=ax.get_yaxis_transform(),
                    ha="center", va="bottom", fontsize=8, zorder=7)
            continue
        xl = P["xl"]; yl = ic + sl*xl; side = P.get("side", 1)   # +1: label above the line, -1: below
        p1, p2 = ax.transData.transform([(xl - 0.1, ic + sl*(xl - 0.1)), (xl + 0.1, ic + sl*(xl + 0.1))])
        ang = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])
        ax.annotate(f"slope = {sgn(sl)}", xy=(xl, yl),
                    xytext=(-4*side*np.sin(ang), 4*side*np.cos(ang)), textcoords="offset points",
                    rotation=np.degrees(ang), rotation_mode="anchor", ha="center",
                    va="bottom" if side > 0 else "top",
                    fontsize=8, zorder=7)
    # solid frame outlines (matplotlib's default "no outline" is written as a transparent state)
    fig.patch.set_edgecolor(fig.patch.get_facecolor())
    for a in fig.axes: a.patch.set_edgecolor(a.patch.get_facecolor())   # panels and color bar
    W, H = fig.get_size_inches(); ax0 = axs[0].get_position()
    print(f"{out}: axes {ax0.width*W:.3f} x {ax0.height*H:.3f} in; "
          f"{3.35/(ax0.width*W):.3f} dex/in, {abs(np.diff(axs[0].get_ylim())[0])/(ax0.height*H):.3f} mag/in")
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
    draw([dict(lz=lz, v=m - 2.5*np.log10(1 - (1 + z)**-2) - 55.3, isQ=isQ, xl=-1.90,
               title=r"Anti-de Sitter ($C_{\rm AdS} = -55.3$)"),
          dict(lz=lz, v=m - 5*lz - 43.0, isQ=isQ, xl=-1.30, title=r"Hubble ($C_{\rm H} = -43.0$)")],
         "Absolute magnitude", "desi_bright_hubble_vs_ads")
    if os.path.exists(sdss):
        zs, ms, isQs = load_sdss(sdss)
        # common magnitude range for both panels; figure size chosen so the panels keep the main
        # figure's scale (1.40 dex/in, 6.80 mag/in), so slopes look equally steep in both figures
        draw([dict(lz=np.log10(zs), v=ms, isQ=isQs, xl=XL_APP[0], side=1, title="SDSS"),
              dict(lz=lz, v=m, isQ=isQ, xl=XL_APP[1], side=1, title="DESI")],
             "Apparent magnitude", "apparent_magnitude_sdss_desi", invert=False, ylim=(11.0, 26.0),
             xarr=XARR_APP, figsize=(5.818, 2.85))
    else:
        print(f"{sdss} not found: appendix figure skipped (run sdss_legacy_selection.py first)")

XL_APP = (("centre", 22.35), ("centre", 22.35))   # slope labels: centred, just above the top band edges
# appendix arrows: opposing arrows share one vertical line; slightly wider spacing than the main
# figure because the DESI band edges bunch up at the faint limit
XARR_APP = {"99.7": -2.45, "95": -2.08, "68": -1.71}

if __name__ == "__main__":
    main(*sys.argv[1:3])
