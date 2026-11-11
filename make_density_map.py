"""
Density-map figure of the paper: the distribution of absolute magnitude at each redshift for
the de Sitter and Hubble (q0 = 1) relations, one panel per sample: (a) CfA + Hewitt & Burbidge
(to 1995), (b) SDSS-I/II (2000-2008), (c) DESI DR1 bright time (2021-2022).
Output: density_map_cfa_sdss_desi.pdf/.png

Each relation is plotted relative to its overall median (all objects of the sample), so both
overall medians lie on the dotted zero line. Color depth is the density of objects in magnitude
at each log10 z, scaled to its maximum in each column of 0.01 in log10 z: green for de Sitter,
reddish purple for Hubble, gray where both are dense. Densities and percentiles use a Gaussian
window in log10 z with standard deviation 0.02, widened in steps up to 0.15 where data are
sparse, until the window holds at least 400 objects (Gaussian-weighted count). In magnitude
they use a Gaussian of 0.15 mag, widened by (2000/N)^(1/4) where that count N is below 2000.
Objects outside the plotted range of log10 z enter the windows near its edges. Solid lines are
running medians, dashed lines the 16th and 84th percentiles. N gives the number of objects within
+-0.125 of each labeled log10 z. The header gives the range spanned by each running median over
the plotted range. All colors are opaque (the overlap color is computed), so the PDF contains no
transparency; the density image is embedded at 400 dpi, and text and lines are vector graphics.

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
from scipy.ndimage import gaussian_filter1d
from make_figure import load_cfa_hb, load_sdss, load_desi, C_ADS, C_H, XLIM, neg   # also sets the house style

COL = {"ds": "#009E73", "h": "#CC79A7"}      # Okabe-Ito bluish green and reddish purple
EDGE = {"ds": "#00664a", "h": "#9b4a7e"}     # darker shades of the same hues, for lines and labels
INK = 0.9                                    # color depth at the maximum of a column
YLIM = (-10.0, 12.0)                         # M - median(M) [mag]; magnitudes increase upward
XB = np.arange(-3.0, 1.0 + 1e-9, 0.01)       # log10 z column edges, wider than XLIM for the edge windows
YB = np.arange(-30.0, 30.0 + 1e-9, 0.05)     # magnitude bin edges [mag]
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


def header(ax, letter, name, years, spans):
    ax.set_title(rf"$\mathbf{{{name}}}$ ({years})", loc="left", pad=15)
    ax.text(-0.075, 1.13, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    tp = lambda c: dict(color=c, fontsize=8)
    box = HPacker(children=[TextArea("medians span ", textprops=tp("0.35")),
                            TextArea(f"de Sitter {spans['ds']:.1f} mag", textprops=tp(EDGE["ds"])),
                            TextArea(", ", textprops=tp("0.35")),
                            TextArea(f"Hubble {spans['h']:.1f} mag", textprops=tp(EDGE["h"]))],
                  align="baseline", pad=0, sep=0)
    ax.add_artist(AnchoredOffsetbox(loc="lower right", child=box, pad=0, borderpad=0, frameon=False,
                                    bbox_to_anchor=(1.0, 1.115), bbox_transform=ax.transAxes))


def main(cfa="cfa_hb_sample.csv.gz", sdss="sdss_legacy_sample.csv.gz", desi="desi_dr1_main_bright_sample.csv.gz",
         out="density_map_cfa_sdss_desi"):
    ads = lambda z, m: m - 2.5*np.log10(1 - (1 + z)**-2) + C_ADS
    hub = lambda z, m: m - 5*np.log10(z) + C_H
    surveys = [("CfA", "to 1995", load_cfa_hb(cfa)), ("SDSS", "2000\u20132008", load_sdss(sdss)),
               ("DESI", "2021\u20132022", load_desi(desi))]
    yc = 0.5*(YB[1:] + YB[:-1]); rv = (yc >= YLIM[0] - 1) & (yc <= YLIM[1] + 1)
    fig, axes = plt.subplots(3, 1, figsize=(5.8, 6.9), sharex=True, gridspec_kw=dict(hspace=0.62))
    print("sample  N  median M_dS  median M_H  span dS  span H  min window count in view")
    for ax, (name, years, (z, m, q)), letter in zip(axes, surveys, "abc"):
        lz = np.log10(z); M = {"ds": ads(z, m), "h": hub(z, m)}
        med = {f: float(np.median(v)) for f, v in M.items()}
        maps = {f: cond_map(lz, v - med[f]) for f, v in M.items()}
        xc = maps["ds"][0]; cv = (xc >= XLIM[0]) & (xc <= XLIM[1])
        img = blend(maps["ds"][1][cv][:, rv], maps["h"][1][cv][:, rv])
        ext = [xc[cv][0] - 0.005, xc[cv][-1] + 0.005, yc[rv][0] - 0.025, yc[rv][-1] + 0.025]
        ax.imshow(img.transpose(1, 0, 2), extent=ext, origin="lower", aspect="auto", interpolation="bilinear", zorder=1)
        spans = {}
        for f in ("h", "ds"):
            qq = maps[f][2]
            ax.plot(xc[cv], qq[0.5][cv], color=EDGE[f], lw=1.3, zorder=5)
            for p in (0.16, 0.84):
                ax.plot(xc[cv], qq[p][cv], color=EDGE[f], lw=0.55, ls=(0, (3, 1.5)), zorder=5)
            spans[f] = float(np.nanmax(qq[0.5][cv]) - np.nanmin(qq[0.5][cv]))
        ax.axhline(0, color="0.15", lw=0.7, ls=(0, (1, 1.6)), zorder=4)
        for x in NPOS:
            ax.text(x, 1.015, kfmt(int(((lz >= x - 0.125) & (lz < x + 0.125)).sum())), transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=7, color="0.4")
        ax.text(-0.012, 1.015, "N", transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="0.4")
        ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
        ax.yaxis.set_major_locator(FixedLocator([-8, -4, 0, 4, 8])); ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: neg(f"{v:+.0f}") if v else "0"))
        ax.xaxis.set_major_locator(MultipleLocator(0.5)); ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.set_ylabel(r"$M - M_{\rm median}$ (mag)")
        header(ax, letter, name, years, spans)
        print(f"{name:5s} {len(z):8d}  {med['ds']:7.2f}  {med['h']:7.2f}  {spans['ds']:5.2f}  {spans['h']:5.2f}  "
              f"{int(min(maps['ds'][3][cv].min(), maps['h'][3][cv].min()))}")
    axes[-1].set_xlabel(r"$\log_{10} z$")
    sw = {k: tuple(blend(np.array([[a]]), np.array([[b]]))[0, 0]) for k, (a, b) in
          {"ds": (1.0, 0.0), "h": (0.0, 1.0), "both": (1.0, 1.0)}.items()}
    handles = [Patch(facecolor=sw["ds"], edgecolor=sw["ds"], lw=0.0, label="de Sitter"),
               Patch(facecolor=sw["h"], edgecolor=sw["h"], lw=0.0, label="Hubble"),
               Patch(facecolor=sw["both"], edgecolor=sw["both"], lw=0.0, label="both"),
               Line2D([0], [0], color="0.15", lw=1.3, label="running median"),
               Line2D([0], [0], color="0.15", lw=0.55, ls=(0, (3, 1.5)), label="16th and 84th percentiles"),
               Line2D([0], [0], color="0.15", lw=0.7, ls=(0, (1, 1.6)), label="overall median")]
    axes[-1].legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.30),
                    handlelength=1.6, columnspacing=1.4)
    # solid frame outlines (matplotlib's default "no outline" is written as a transparent state), as in make_figure.py
    fig.patch.set_edgecolor(fig.patch.get_facecolor())
    for a in fig.axes:
        a.patch.set_edgecolor(a.patch.get_facecolor())
    fig.savefig(out + ".pdf", dpi=400)
    fig.savefig(out + ".png", dpi=300)


if __name__ == "__main__":
    main()
