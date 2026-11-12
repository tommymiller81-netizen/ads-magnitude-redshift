"""
Density-map figures of the paper, one panel per sample: (a) CfA + Hewitt & Burbidge (to 1995),
(b) SDSS-I/II (2000-2008), (c) DESI DR1 bright time (2021-2022).

1. density_map_cfa_sdss_desi.pdf/.png: the distribution of absolute magnitude at each redshift for
   the de Sitter and Hubble (q0 = 1) relations, each plotted relative to its overall median (all
   objects of the sample), so both overall medians lie on the dotted zero line. Green is de Sitter,
   reddish purple Hubble, gray where both are dense.
2. apparent_density_map_cfa_sdss_desi.pdf/.png (appendix): the same for the apparent magnitude m.

The script also prints the statistics of the paper's table: for M (de Sitter, Hubble) and m, the
overall median, the range spanned by the running median over the plotted range, and the unweighted
least-squares slope against log10 z (all objects); for DESI also the Hubble values for q0 = 0.5
(Mattig luminosity distance) and flat LambdaCDM (Omega_m = 0.3). Values are quoted to one decimal.

Method. Color depth is the density of objects in magnitude at each log10 z, scaled to its maximum
in each column of 0.01 in log10 z. Densities and percentiles use a Gaussian window in log10 z with
standard deviation 0.02, widened in steps up to 0.15 where data are sparse, until the window holds
at least 400 objects (weighted by the window). In magnitude they use a Gaussian of 0.15 mag, widened
by (2000/N)^(1/4) where that weighted count N is below 2000. Objects outside the plotted range of
log10 z enter the windows near its edges. Solid lines are running medians, dashed lines the 16th and
84th percentiles. N gives the number of objects within +-0.125 of each labeled log10 z. A short tick
marks the log10 z above which quasars are the majority (quasar fraction in bins of 0.02 in log10 z,
smoothed with a Gaussian of 0.05). All colors are opaque (the overlap color is computed), so the PDFs
contain no transparency; the density images are embedded at 400 dpi, text and lines are vector.

Usage:  python make_density_map.py      (inputs as for make_figure.py, in the working directory)
Requires: numpy, pandas, scipy, matplotlib, and make_figure.py (loaders, constants, house style).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # make_figure.py sits next to this script
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, MultipleLocator, FuncFormatter
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
from make_figure import load_cfa_hb, load_sdss, load_desi, C_ADS, C_H, XLIM, neg   # also sets the house style

COL = {"ds": "#009E73", "h": "#CC79A7", "m": "#5f5f5f"}   # Okabe-Ito bluish green and reddish purple; gray for m
EDGE = {"ds": "#00664a", "h": "#9b4a7e", "m": "#1a1a1a"}  # darker shades of the same hues, for lines and labels
INK = 0.9                                    # color depth at the maximum of a column
YLIM = (-10.0, 12.0)                         # figure 1: M - median(M) [mag]; magnitudes increase upward
YLIM_M = (10.0, 25.0)                        # figure 2: m [mag]
XB = np.arange(-3.0, 1.0 + 1e-9, 0.01)       # log10 z column edges, wider than XLIM for the edge windows
YB = np.arange(-30.0, 30.0 + 1e-9, 0.05)     # magnitude bin edges [mag], relative to the overall median
SIGS = (0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.125, 0.15)   # window widths tried, in log10 z
NMIN, SY0, NREF = 400, 0.15, 2000            # min. objects per window; magnitude smoothing [mag] and its reference count
NPOS = np.arange(-2.5, 0.5 + 1e-9, 0.25)     # positions of the N labels (objects within +-0.125)


def kfmt(n):
    """Compact count: 514, 1.3k, 33k, 2.6M."""
    return f"{n/1e6:.1f}M" if n >= 1e6 else f"{n/1e3:.0f}k" if n >= 1e4 else f"{n/1e3:.1f}k" if n >= 1e3 else str(n)


def cond_map(lz, dm):
    """Density of dm in each log10 z column (scaled to the column maximum), running percentiles,
    and the Gaussian-weighted count of each window."""
    xc = 0.5*(XB[1:] + XB[:-1])
    H = np.histogram2d(lz, dm, bins=[XB, YB])[0]; n = H.sum(1)
    W = np.empty((len(xc), len(xc))); neff = np.empty(len(xc))
    for i, x0 in enumerate(xc):
        for s in SIGS:
            w = np.exp(-0.5*((xc - x0)/s)**2)
            if (w*n).sum() >= NMIN:
                break
        W[i], neff[i] = w, (w*n).sum()
    S = W @ H
    for i in range(len(xc)):
        sy = SY0*max(1.0, (NREF/max(neff[i], 1.0))**0.25)
        S[i] = gaussian_filter1d(S[i], sy/(YB[1] - YB[0]), mode="constant")
    tot = S.sum(1); ok = tot > 0
    cdf = np.concatenate([np.zeros((len(xc), 1)), np.cumsum(S, 1)], 1)/np.where(ok, tot, 1.0)[:, None]
    ramp = 1e-12*np.arange(len(YB))                              # makes the cumulative sum strictly increasing
    q = {p: np.array([np.interp(p, cdf[i] + ramp, YB) if ok[i] else np.nan for i in range(len(xc))])
         for p in (0.16, 0.5, 0.84)}
    D = np.where(ok[:, None], S/np.maximum(S.max(1, keepdims=True), 1e-300), 0.0)
    return xc, D, q, neff


def blend(Dd, Dh):
    """Opaque two-color image: hue is the density-weighted mean of the two colors, depth is
    1 - (1 - Dd)(1 - Dh). Equal densities give the mean hue (gray), so neither relation dominates."""
    cd, ch = np.array(to_rgb(COL["ds"])), np.array(to_rgb(COL["h"]))
    Dd, Dh = np.clip(Dd, 0, 1), np.clip(Dh, 0, 1)
    hue = (Dd[..., None]*cd + Dh[..., None]*ch)/(Dd + Dh + 1e-12)[..., None]
    depth = (1 - (1 - Dd)*(1 - Dh))[..., None]
    return 1 - INK*depth*(1 - hue)


def blend1(D, color):
    """Opaque one-color image."""
    return 1 - INK*np.clip(D, 0, 1)[..., None]*(1 - np.array(to_rgb(color)))


def qso_majority(lz, isq):
    """log10 z above which quasars are the majority: quasar fraction in bins of 0.02, smoothed with
    a Gaussian of 0.05; start of the last run of bins with fraction >= 1/2 (None if there is none)."""
    e = np.arange(XLIM[0], XLIM[1] + 0.01, 0.02); c = 0.5*(e[1:] + e[:-1])
    i = np.digitize(lz, e) - 1; ok = (i >= 0) & (i < len(c))
    n = gaussian_filter1d(np.bincount(i[ok], minlength=len(c)).astype(float), 2.5)
    nq = gaussian_filter1d(np.bincount(i[ok], weights=isq[ok].astype(float), minlength=len(c)), 2.5)
    above = nq >= 0.5*n
    if not above[-1]:
        return None
    k = len(c) - np.argmax(~above[::-1])
    return float(c[k]) if k < len(c) else None


def mark_quasars(ax, x50, where):
    """Short tick at x50 with 'galaxies' to its left and 'quasars' to its right, at the top or bottom."""
    tr = ax.get_xaxis_transform()
    y0, y1, yt = (0.88, 1.0, 0.94) if where == "top" else (0.0, 0.12, 0.06)
    ax.plot([x50, x50], [y0, y1], transform=tr, color="0.3", lw=0.7, zorder=6, solid_capstyle="butt")
    kw = dict(transform=tr, va="center", fontsize=7, color="0.3", zorder=6)
    ax.text(x50 - 0.03, yt, "galaxies", ha="right", **kw)
    ax.text(x50 + 0.03, yt, "quasars", ha="left", **kw)


def header(ax, letter, name, years, parts):
    """Panel letter, sample name and years on the left; a colored statistics line on the right."""
    ax.set_title(rf"$\mathbf{{{name}}}$ ({years})", loc="left", pad=15)
    ax.text(-0.075, 1.13, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    box = HPacker(children=[TextArea(t, textprops=dict(color=c, fontsize=8)) for t, c in parts],
                  align="baseline", pad=0, sep=0)
    ax.add_artist(AnchoredOffsetbox(loc="lower right", child=box, pad=0, borderpad=0, frameon=False,
                                    bbox_to_anchor=(1.0, 1.115), bbox_transform=ax.transAxes))


def frame(ax, lz, ylim, yticks, yfmt, ylabel):
    """N row, axis limits, ticks and labels shared by both figures."""
    for x in NPOS:
        ax.text(x, 1.015, kfmt(int(((lz >= x - 0.125) & (lz < x + 0.125)).sum())), transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=7, color="0.4")
    ax.text(-0.012, 1.015, "N", transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="0.4")
    ax.set_xlim(*XLIM); ax.set_ylim(*ylim)
    ax.yaxis.set_major_locator(yticks); ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_formatter(FuncFormatter(yfmt))
    ax.xaxis.set_major_locator(MultipleLocator(0.5)); ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.set_ylabel(ylabel)


def image(ax, xc, rgb_cols, yc, rows, yoff=0.0):
    """Place an image of columns xc (within XLIM) and magnitude rows yc (offset by yoff)."""
    ext = [xc[0] - 0.005, xc[-1] + 0.005, yc[rows][0] - 0.025 + yoff, yc[rows][-1] + 0.025 + yoff]
    ax.imshow(rgb_cols.transpose(1, 0, 2), extent=ext, origin="lower", aspect="auto", interpolation="bilinear", zorder=1)


def save(fig, out):
    # solid frame outlines (matplotlib's default "no outline" is written as a transparent state), as in make_figure.py
    fig.patch.set_edgecolor(fig.patch.get_facecolor())
    for a in fig.axes:
        a.patch.set_edgecolor(a.patch.get_facecolor())
    fig.savefig(out + ".pdf", dpi=400)
    fig.savefig(out + ".png", dpi=300)


def main(cfa="cfa_hb_sample.csv.gz", sdss="sdss_legacy_sample.csv.gz", desi="desi_dr1_main_bright_sample.csv.gz"):
    ads = lambda z, m: m - 2.5*np.log10(1 - (1 + z)**-2) + C_ADS
    hub = lambda z, m: m - 5*np.log10(z) + C_H
    zg = np.linspace(0, 8, 80001)                                   # flat LambdaCDM, Omega_m = 0.3
    chi = cumulative_trapezoid(1/np.sqrt(0.3*(1 + zg)**3 + 0.7), zg, initial=0)
    dl_alt = {"q0 = 0.5": lambda z: 2*(1 + z - np.sqrt(1 + z)),     # Mattig, D_L H0/c
              "LCDM, Omega_m = 0.3": lambda z: (1 + z)*np.interp(z, zg, chi)}
    surveys = [("CfA", "to 1995", load_cfa_hb(cfa)), ("SDSS", "2000\u20132008", load_sdss(sdss)),
               ("DESI", "2021\u20132022", load_desi(desi))]
    yc = 0.5*(YB[1:] + YB[:-1])
    kw = dict(figsize=(5.8, 6.9), sharex=True, gridspec_kw=dict(hspace=0.62))
    fig1, axes1 = plt.subplots(3, 1, **kw); fig2, axes2 = plt.subplots(3, 1, **kw)
    stats = []
    for ax1, ax2, (name, years, (z, m, isq)), letter in zip(axes1, axes2, surveys, "abc"):
        lz = np.log10(z); x50 = qso_majority(lz, isq)
        V = {"ds": ads(z, m), "h": hub(z, m), "m": m}
        med = {f: float(np.median(v)) for f, v in V.items()}
        maps = {f: cond_map(lz, v - med[f]) for f, v in V.items()}
        xc = maps["ds"][0]; cv = (xc >= XLIM[0]) & (xc <= XLIM[1])
        span = {f: float(np.nanmax(maps[f][2][0.5][cv]) - np.nanmin(maps[f][2][0.5][cv])) for f in V}
        for f in V:
            stats.append((name, f, len(z), isq.mean(), med[f], span[f], np.polyfit(lz, V[f], 1)[0], x50))
        if name == "DESI":
            for lab, dl in dl_alt.items():
                v = m - 5*np.log10(dl(z)) + C_H; mp = cond_map(lz, v - np.median(v))
                stats.append((name, f"h, {lab}", len(z), isq.mean(), float(np.median(v)),
                              float(np.nanmax(mp[2][0.5][cv]) - np.nanmin(mp[2][0.5][cv])), np.polyfit(lz, v, 1)[0], x50))
        # figure 1: both relations relative to their overall medians
        rv = (yc >= YLIM[0] - 1) & (yc <= YLIM[1] + 1)
        image(ax1, xc[cv], blend(maps["ds"][1][cv][:, rv], maps["h"][1][cv][:, rv]), yc, rv)
        for f in ("h", "ds"):
            qq = maps[f][2]
            ax1.plot(xc[cv], qq[0.5][cv], color=EDGE[f], lw=1.3, zorder=5)
            for p in (0.16, 0.84):
                ax1.plot(xc[cv], qq[p][cv], color=EDGE[f], lw=0.55, ls=(0, (3, 1.5)), zorder=5)
        ax1.axhline(0, color="0.15", lw=0.7, ls=(0, (1, 1.6)), zorder=4)
        frame(ax1, lz, YLIM, FixedLocator([-8, -4, 0, 4, 8]), lambda v, _: neg(f"{v:+.0f}") if v else "0",
              r"$M - M_{\rm median}$ (mag)")
        header(ax1, letter, name, years, [("medians span ", "0.35"), (f"de Sitter {span['ds']:.1f} mag", EDGE["ds"]),
                                          (", ", "0.35"), (f"Hubble {span['h']:.1f} mag", EDGE["h"])])
        if x50 is not None:
            mark_quasars(ax1, x50, "top")
        # figure 2: apparent magnitude
        rv2 = (yc + med["m"] >= YLIM_M[0] - 1) & (yc + med["m"] <= YLIM_M[1] + 1)
        image(ax2, xc[cv], blend1(maps["m"][1][cv][:, rv2], COL["m"]), yc, rv2, yoff=med["m"])
        qq = maps["m"][2]
        ax2.plot(xc[cv], qq[0.5][cv] + med["m"], color=EDGE["m"], lw=1.3, zorder=5)
        for p in (0.16, 0.84):
            ax2.plot(xc[cv], qq[p][cv] + med["m"], color=EDGE["m"], lw=0.55, ls=(0, (3, 1.5)), zorder=5)
        frame(ax2, lz, YLIM_M, MultipleLocator(2), lambda v, _: f"{v:.0f}", r"apparent magnitude $m$")
        header(ax2, letter, name, years, [("median spans ", "0.35"), (f"{span['m']:.1f} mag", "0.15")])
        if x50 is not None:
            mark_quasars(ax2, x50, "bottom")
    for axes in (axes1, axes2):
        axes[-1].set_xlabel(r"$\log_{10} z$")
    sw = {k: tuple(blend(np.array([[a]]), np.array([[b]]))[0, 0]) for k, (a, b) in
          {"ds": (1.0, 0.0), "h": (0.0, 1.0), "both": (1.0, 1.0)}.items()}
    lines = [Line2D([0], [0], color="0.15", lw=1.3, label="running median"),
             Line2D([0], [0], color="0.15", lw=0.55, ls=(0, (3, 1.5)), label="16th and 84th percentiles")]
    h1 = [Patch(facecolor=sw[k], edgecolor=sw[k], lw=0.0, label=lab) for k, lab in
          (("ds", "de Sitter"), ("h", "Hubble"), ("both", "both"))] + lines + \
         [Line2D([0], [0], color="0.15", lw=0.7, ls=(0, (1, 1.6)), label="overall median")]
    sm = tuple(blend1(np.array([[1.0]]), COL["m"])[0, 0])
    h2 = [Patch(facecolor=sm, edgecolor=sm, lw=0.0, label="density")] + lines
    axes1[-1].legend(handles=h1, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.30), handlelength=1.6, columnspacing=1.4)
    axes2[-1].legend(handles=h2, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.30), handlelength=1.6, columnspacing=1.4)
    save(fig1, "density_map_cfa_sdss_desi"); save(fig2, "apparent_density_map_cfa_sdss_desi")
    names = {"ds": "de Sitter", "h": "Hubble (q0 = 1)", "m": "apparent m"}
    print(f"{'sample':6s} {'quantity':30s} {'N':>9s} {'QSO %':>6s} {'median':>7s} {'span':>5s} {'slope':>6s}  quasar majority above log10 z")
    for name, f, n, fq, md, sp, sl, x50 in stats:
        lab = names.get(f, "Hubble (" + f[3:] + ")")
        print(f"{name:6s} {lab:30s} {n:9d} {100*fq:6.0f} {md:7.1f} {sp:5.1f} {neg(f'{sl:+.1f}'):>6s}  {x50:+.2f}")


if __name__ == "__main__":
    main()
