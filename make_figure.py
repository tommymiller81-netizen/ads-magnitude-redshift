"""
Reproduce the figures of the paper.

1. Main figure: absolute magnitude against log10 z. Rows: (a) CfA + Hewitt & Burbidge (to 1995),
   (b) SDSS-I/II (2000-2008), (c) DESI DR1 bright time (2021-2022).
   Columns: (1) anti-de Sitter, (2) Hubble (q0 = 1).
   Output: absolute_magnitude_cfa_sdss_desi.pdf/.png
2. Appendix figure: apparent magnitude m against log10 z, (a) CfA, (b) SDSS, (c) DESI.
   Output: apparent_magnitude_cfa_sdss_desi.pdf/.png
Inputs: cfa_hb_sample.csv.gz (cfa_hb_selection.py), sdss_legacy_sample.csv.gz
        (sdss_legacy_selection.py) and desi_dr1_main_bright_sample.csv.gz (desi_bright_selection.py).

Panel letters name the sample (a = CfA, b = SDSS, c = DESI) in both figures; digits name the
relation. Magnitudes increase upward, so fainter objects plot higher. Gray shading is the
fraction of each sample per cell, on one logarithmic scale for all panels, so samples of different
size share one scale; cells holding a single object are drawn in one fixed light gray. Each panel
shows 18 mag; all panels in each figure have the same physical size (1.50 x 1.665 in; 2.23 dex/in,
10.8 mag/in), so slopes can be compared by eye across all panels. Labels give the least-squares
slope (mag per unit log10 z) and the full widths of the central 68%, 95% and 99.7% of the sample.
Use the PDFs for JHEP (no transparency); the PNGs are previews. The figure captions are in the
paper (LaTeX), not in the images.

Usage:  python make_figure.py
Requires: numpy, pandas, matplotlib.
"""
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm

mpl.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.titlelocation": "left", "lines.linewidth": 1.2, "patch.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42})

C_ADS, C_H = -55.3, -43.0                   # constants of the anti-de Sitter and Hubble relations (r band)
XLIM = (-2.6, 0.75)                         # log10 z range of every panel
YMAG = 18.0                                 # magnitude range of every panel
PANEL_W = 1.50                              # panel width [in]; height = PANEL_W*1.11025 (4.84 mag per dex on paper)
NMIN_G, NMIN_Q = 100, 100                   # min. galaxies / quasars per 0.1-dex bin for running medians
MINCNT = 1                                  # cells with at least one object are shaded
SINGLE_GRAY = (0.88, 0.88, 0.88, 1.0)       # cells with a single object
BANDC = ["#fffbe6", "#fff3bf", "#ffe88a"]; EDGEC = "#dfc77d"        # 99.7%, 95%, 68% bands; opaque
CG, CQ, CF = "#0072B2", "#D55E00", "#000000"                         # Okabe-Ito blue / vermillion; fit
FIT_DASH = tuple(v/1.1 for v in (4.0, 2.0))  # least-squares line: 4-pt dashes, 2-pt gaps (line width 1.1 pt)
SLOPE_FS, SLOPE_D, WGAP, LABEL_X = 9, 2.4, 1.25, 0.97   # slope label (pt; mag below top edge); widths line; x
neg = lambda s: s.replace("-", "\u2212")                 # typographic minus
sgn = lambda x: neg(f"{x:+.1f}")                         # one decimal, explicit sign

def panel(ax, P, edges, cen):
    """One panel: percentile bands, weighted hexbin, least-squares line, running medians.
    P: dict with lz (log10 z), v (magnitude), isQ (quasar flag), ylim, title (or None).
    Returns (hexbin, slope)."""
    lz, v, isQ = P["lz"], P["v"], P["isQ"]; idx = np.digitize(lz, edges) - 1
    p = np.percentile(v, [0.15, 2.5, 16, 50, 84, 97.5, 99.85])
    for (lo, hi), c in zip([(0, 6), (1, 5), (2, 4)], BANDC):
        ax.axhspan(p[lo], p[hi], color=c, lw=0, zorder=0)
    for y in (p[0], p[1], p[2], p[4], p[5], p[6]):
        ax.axhline(y, color=EDGEC, lw=0.4, zorder=1)
    y0, y1 = P["ylim"]
    ny = int(round(70*(y1 - y0)/YMAG))
    ext = (*XLIM, y0, y1)
    w = np.full(len(v), 1.0/len(v))                       # fraction of the sample per cell
    hb = ax.hexbin(lz, v, C=w, reduce_C_function=np.sum, gridsize=(110, ny), extent=ext, norm=LogNorm(),
                   cmap="Greys", mincnt=MINCNT, linewidths=0, zorder=2)
    hc = ax.hexbin(lz, v, gridsize=(110, ny), extent=ext, mincnt=MINCNT, visible=False)
    hb._counts = np.asarray(hc.get_array()); hc.remove()  # objects per cell, same cell order
    assert len(hb._counts) == len(hb.get_array())
    sl, ic = np.polyfit(lz, v, 1); xx = np.array(XLIM)
    ax.plot(xx, ic + sl*xx, color=CF, lw=1.1, ls=(0, FIT_DASH), zorder=4)
    def bmed(sel, nmin):
        vv, ii = v[sel], idx[sel]
        return np.array([np.median(vv[ii == b]) if (ii == b).sum() >= nmin else np.nan for b in range(len(cen))])
    medG, medQ = bmed(~isQ, NMIN_G), bmed(isQ, NMIN_Q)
    ax.plot(cen, medQ, "-", color=CQ, lw=0.9, zorder=9)
    ax.plot(cen, medQ, ls="none", marker="s", mfc=CQ, mec=CQ, mew=0.3, ms=3.0, zorder=9)
    ax.plot(cen, medG, "o-", color=CG, ms=2.6, lw=1.1, zorder=10)
    ax.set_ylim(y0, y1); ax.set_xlim(*XLIM); ax.yaxis.set_major_locator(MultipleLocator(2))
    if P.get("title"): ax.set_title(P["title"])
    # labels: slope (headline) and the band widths, right-aligned in the upper right
    ax.text(LABEL_X, y1 - SLOPE_D, f"slope = {sgn(sl)}", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=SLOPE_FS, zorder=7)
    ax.text(LABEL_X, y1 - SLOPE_D - WGAP, f"widths {p[4]-p[2]:.1f}/{p[5]-p[1]:.1f}/{p[6]-p[0]:.1f} mag",
            transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=7, color="0.35", zorder=7)
    return hb, sl

def draw_grid(rows, out, ylabels, letters, hide_inner_y=True, legend_ncol=3):
    """Panels in a grid: rows[r][c] is a panel dict; ylabels[r] labels row r (first panel);
    letters[r][c] is the panel letter. Legend and color bar share one strip below the panels."""
    edges = np.arange(XLIM[0], XLIM[1], 0.1); cen = 0.5*(edges[1:] + edges[:-1])
    nr, nc = len(rows), len(rows[0]); target = (PANEL_W, PANEL_W*1.11025)
    fig, axs = plt.subplots(nr, nc, figsize=(5.85, 3.0*nr), constrained_layout=True, squeeze=False)
    hbs = []
    for r, row in enumerate(rows):
        for c, P in enumerate(row):
            ax = axs[r, c]
            hb, sl = panel(ax, dict(P), edges, cen); hbs.append(hb)
            if r == nr - 1: ax.set_xlabel(r"$\log_{10} z$")
            else: ax.tick_params(labelbottom=False)
            if c == 0: ax.set_ylabel(ylabels[r])
            elif hide_inner_y: ax.tick_params(labelleft=False)
            ax.text(-0.18, 1.02, letters[r][c], transform=ax.transAxes, fontweight="bold", fontsize=10,
                    va="bottom", ha="left")
    h = [Patch(facecolor=BANDC[2], edgecolor=EDGEC, lw=0.4, label="central 68%"),
         Line2D([], [], color=CG, marker="o", ms=2.6, lw=1.1, label="galaxy median"),
         Patch(facecolor=BANDC[1], edgecolor=EDGEC, lw=0.4, label="central 95%"),
         Line2D([], [], color=CQ, marker="s", ms=3.0, mfc=CQ, mec=CQ, mew=0.3, lw=0.9, label="quasar median"),
         Patch(facecolor=BANDC[0], edgecolor=EDGEC, lw=0.4, label="central 99.7%"),
         Line2D([], [], color=CF, lw=1.1, ls=(0, FIT_DASH), label="least-squares fit")]
    if legend_ncol == 2: h = [h[0], h[2], h[4], h[1], h[3], h[5]]   # column 1 bands, column 2 lines
    leg = fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.0, 0.0), ncol=legend_ncol, fontsize=8,
                     columnspacing=0.9, handletextpad=0.5, borderaxespad=0.0)
    # one logarithmic gray scale for all panels; single-object cells in a fixed light gray (otherwise a
    # single object in a small sample is drawn as dark as hundreds in a large one)
    lo = min(h_.get_clim()[0] for h_ in hbs); hi = max(h_.get_clim()[1] for h_ in hbs)
    for h_ in hbs:
        col = plt.get_cmap("Greys")(LogNorm(lo, hi)(np.asarray(h_.get_array())))
        col[h_._counts == 1] = SINGLE_GRAY
        h_.set_array(None); h_.set_facecolor(col)
    sm = cm.ScalarMappable(norm=LogNorm(lo, hi), cmap="Greys")
    # panels at the common printed size: adjust the figure size until the axes match the target
    for _ in range(6):
        fig.canvas.draw(); b_ = axs[0, 0].get_position(); W, H = fig.get_size_inches()
        fig.set_size_inches(W + nc*(target[0] - b_.width*W), H + nr*(target[1] - b_.height*H))
    # legend and color bar side by side; the color bar label ends at the right edge of the panels
    fig.canvas.draw(); W, H = fig.get_size_inches()
    x_ax1 = axs[0, nc - 1].get_position().x1*W
    gap, cbw, cbh = 0.25, 1.15, 0.086
    bb = leg.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    cy = bb.y0 + 0.62*bb.height
    cax = fig.add_axes([0.5, cy/H, cbw/W, cbh/H])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("fraction of sample per cell", fontsize=8)
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
    print(f"{out}: axes {b_.width*W:.3f} x {b_.height*H:.3f} in; {(XLIM[1]-XLIM[0])/(b_.width*W):.3f} dex/in, "
          f"{YMAG/(b_.height*H):.3f} mag/in; color scale {lo:.2e} to {hi:.2e}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}")
    plt.close(fig); print(f"wrote {out}.png/.pdf")

def load_cfa_hb(path):
    c = pd.read_csv(path)
    return c.z.to_numpy(float), c.Bmag.to_numpy(float), (c["sample"] == "QSO").to_numpy()

def load_sdss(path):
    s = pd.read_csv(path)
    return s.z_use.to_numpy(float), s.modelMag_r.to_numpy(float), (s["sample"] == "QSO").to_numpy()

def load_desi(path):
    d = pd.read_csv(path)
    z = d.z.to_numpy(float); m = 22.5 - 2.5*np.log10(d.flux_r.to_numpy(float))
    return z, m, (d.spectype == "QSO").to_numpy()

def main(cfa="cfa_hb_sample.csv.gz", sdss="sdss_legacy_sample.csv.gz", desi="desi_dr1_main_bright_sample.csv.gz"):
    yl_H = (-29.5, -11.5); yl_A = (yl_H[0] - 14.0, yl_H[1] - 14.0)   # 18 mag each
    ads = lambda z, m: m - 2.5*np.log10(1 - (1 + z)**-2) + C_ADS
    hub = lambda z, m: m - 5*np.log10(z) + C_H
    # observing periods: CfA June 1995 catalog and Hewitt & Burbidge (1993), SDSS-I/II spectra, DESI DR1
    surveys = [("CfA", "to 1995", load_cfa_hb(cfa)), ("SDSS", "2000\u20132008", load_sdss(sdss)),
               ("DESI", "2021\u20132022", load_desi(desi))]
    rows, labels, app = [], [], []
    for k, (name, years, (z, m, q)) in enumerate(surveys):
        lz = np.log10(z)
        rows.append([dict(lz=lz, v=ads(z, m), isQ=q, ylim=yl_A, title="Anti-de Sitter" if k == 0 else None),
                     dict(lz=lz, v=hub(z, m), isQ=q, ylim=yl_H, title="Hubble" if k == 0 else None)])
        labels.append(rf"$\mathbf{{{name}}}$ ({years})" + "\nabsolute magnitude $M$")
        app.append(dict(lz=lz, v=m, isQ=q, ylim=(10.0, 28.0), title=f"{name} ({years})"))
    draw_grid(rows, "absolute_magnitude_cfa_sdss_desi", labels, [["a1", "a2"], ["b1", "b2"], ["c1", "c2"]],
              hide_inner_y=False, legend_ncol=2)
    draw_grid([app], "apparent_magnitude_cfa_sdss_desi", ["Apparent magnitude $m$"], [["a", "b", "c"]],
              hide_inner_y=True, legend_ncol=3)

if __name__ == "__main__":
    main()
