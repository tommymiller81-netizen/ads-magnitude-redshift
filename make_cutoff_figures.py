"""
Figures for the redshift-cutoff subsection and appendix A of the paper. No data are needed.

1. cutoff_number_redshift: predicted number of objects against redshift z for a uniform
   population in the static de Sitter metric. dN/dz is proportional to (z+1)^-3 [z(z+2)]^(1/2),
   scaled to a peak of 1; N_z, the fraction of all objects below z, is
   [arctan(sqrt(z/(z+2))) - sqrt(z(z+2)) / (2 (z+1)^2)] / (pi/4).
   Mode (3/2)^(1/2) - 1 = 0.2247, median 1.4754.
2. volume_area_circumference: volume V/R^3, surface area S/R^2 and circumference C/R of a sphere
   centered on the observer, against its radial coordinate in each system: pseudo-Euclidean
   r1/R, elliptical chi, hyperbolic h/R, with Euclidean values for comparison.
   V = 2 pi R^3 (chi - sin chi cos chi), S = 4 pi R^2 sin^2 chi, C = 2 pi R sin chi, where
   sin chi = r1/R = tanh(h/R). Horizon values: pi^2, 4 pi, 2 pi.

Use the PDFs for JHEP; the PNGs are previews. The captions are in the paper (LaTeX).

Usage:  python make_cutoff_figures.py
Requires: numpy, matplotlib (house style from make_figure.py).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import make_figure  # noqa: F401  (sets the house style: fonts, open frame, frameless legends)

BLUE, VERMILION = "#0072B2", "#D55E00"                          # Okabe-Ito, color-blind safe
C_EUC, C_PE, C_ELL, C_HYP = "0.45", "#E69F00", "#009E73", "#0072B2"
DASH = (0, (3, 2))
LABEL_FS = 8


def dndz(z):
    return (z + 1.0)**-3 * np.sqrt(z * (z + 2.0))


def n_below(z):
    return np.arctan(np.sqrt(z / (z + 2.0))) - np.sqrt(z * (z + 2.0)) / (2.0 * (z + 1.0)**2)


Z_MODE = np.sqrt(1.5) - 1.0
N_TOTAL = np.pi / 4.0                                            # n_below(infinity)


def median_redshift():
    lo, hi = 0.0, 10.0
    for _ in range(80):                                          # bisection on N_z = N_total/2
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if n_below(mid) < 0.5 * N_TOTAL else (lo, mid)
    return 0.5 * (lo + hi)


Z_MEDIAN = median_redshift()


def save(fig, out):
    # solid frame outlines (matplotlib's default "no outline" is written as a transparent state)
    fig.patch.set_edgecolor(fig.patch.get_facecolor())
    for a in fig.axes:
        a.patch.set_edgecolor(a.patch.get_facecolor())
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=300)


def tint(color, amount):
    """Opaque mix of a hex color with white (no transparency in the PDF)."""
    c = np.array([int(color[i:i + 2], 16) for i in (1, 3, 5)]) / 255.0
    return tuple(amount * c + (1.0 - amount))


def cutoff_figure(zmax=5.0, integral=False):
    """dN/dz normalized to unit area, shaded in two equal halves at the median.
    integral=True adds a lower panel with N_z, the fraction of all objects below z."""
    z = np.linspace(0.0, zmax, 2501)
    p = dndz(z) / N_TOTAL                                        # unit area over 0 < z < infinity
    p_mode, p_median = dndz(Z_MODE) / N_TOTAL, dndz(Z_MEDIAN) / N_TOTAL
    beyond = 1.0 - n_below(zmax) / N_TOTAL
    if integral:
        fig, (ax, ax2) = plt.subplots(2, 1, figsize=(4.3, 3.6), sharex=True,
                                      gridspec_kw=dict(height_ratios=[2.2, 1.25], hspace=0.12))
    else:
        fig, ax = plt.subplots(figsize=(4.3, 2.4))
    lo, hi = z <= Z_MEDIAN, z >= Z_MEDIAN
    ax.fill_between(z[lo], 0, p[lo], color=tint(BLUE, 0.38), lw=0, zorder=1)
    ax.fill_between(z[hi], 0, p[hi], color=tint(BLUE, 0.14), lw=0, zorder=1)
    ax.plot(z, p, color=BLUE, lw=1.4, zorder=3)
    ax.plot([Z_MEDIAN, Z_MEDIAN], [0, p_median], color=BLUE, lw=0.8, zorder=2)
    ax.plot([Z_MODE], [p_mode], "o", ms=3.2, color=BLUE, mec=BLUE, mew=0.5, zorder=4)
    ax.text(Z_MODE + 0.09, p_mode + 0.004, f"mode {Z_MODE:.2f}", fontsize=LABEL_FS, color=BLUE,
            ha="left", va="bottom")
    ax.text(Z_MEDIAN + 0.07, p_median + 0.012, f"median {Z_MEDIAN:.2f}", fontsize=LABEL_FS,
            color=BLUE, ha="left", va="bottom")
    ax.text(0.72, 0.07, "half", fontsize=LABEL_FS, color="0.15", ha="center", va="center")
    ax.text(2.05, 0.045, "half", fontsize=LABEL_FS, color="0.15", ha="center", va="center")
    ax.annotate(f"{100 * beyond:.0f}% beyond $z$ = {zmax:.0f}", xy=(zmax, 0.025), xytext=(zmax - 0.9, 0.025),
                fontsize=LABEL_FS, color="0.15", ha="right", va="center",
                arrowprops=dict(arrowstyle="->", lw=0.6, color="0.25", shrinkA=2, shrinkB=0))
    ax.set_xlim(0, zmax)
    ax.set_ylim(0, 0.56)
    ax.xaxis.set_major_locator(MultipleLocator(1)); ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.set_ylabel(r"fraction of objects per unit $z$")
    if not integral:
        ax.set_xlabel(r"redshift $z$")
        return fig
    frac = n_below(z) / N_TOTAL
    ax2.plot(z, frac, color=BLUE, lw=1.4, zorder=3)
    ax2.plot([0, Z_MEDIAN, Z_MEDIAN], [0.5, 0.5, 0], color=BLUE, lw=0.7, ls=DASH, zorder=2)
    ax2.text(zmax - 0.03, frac[-1] + 0.07, f"{100 * (1 - beyond):.0f}% below $z$ = {zmax:.0f}", fontsize=LABEL_FS,
             color="0.15", ha="right", va="bottom")
    ax2.set_ylim(0, 1.05)
    ax2.yaxis.set_major_locator(MultipleLocator(0.5)); ax2.yaxis.set_minor_locator(MultipleLocator(0.25))
    ax2.set_ylabel("fraction below $z$")
    ax2.set_xlabel(r"redshift $z$")
    fig.align_ylabels([ax, ax2])
    return fig


def sphere(chi):
    """Volume / R^3, surface area / R^2 and circumference / R inside the angle chi."""
    return 2*np.pi*(chi - np.sin(chi)*np.cos(chi)), 4*np.pi*np.sin(chi)**2, 2*np.pi*np.sin(chi)


def volume_figure(out="volume_area_circumference", xmax=4.0):
    chi = np.linspace(0.0, np.pi/2, 1201)                        # elliptical: x = chi
    r1 = np.sin(chi)                                             # pseudo-Euclidean: x = r1/R
    h = np.linspace(0.0, xmax, 1601)                             # hyperbolic: x = h/R
    chi_h = np.arctan(np.sinh(h))                                # chi = gd(h/R)
    r = np.linspace(0.0, xmax, 801)                              # Euclidean: x = r/R
    euclid = (4*np.pi/3*r**3, 4*np.pi*r**2, 2*np.pi*r)
    panels = [("volume / R$^3$", r"$\pi^2$", np.pi**2, 14.0, 2),
              ("surface area / R$^2$", r"$4\pi$", 4*np.pi, 17.0, 5),
              ("circumference / R", r"$2\pi$", 2*np.pi, 8.5, 2)]
    fig, axes = plt.subplots(1, 3, figsize=(5.8, 2.15))
    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.27, top=0.97, wspace=0.42)
    for k, (ax, (ylab, hlab, hval, ytop, ystep)) in enumerate(zip(axes, panels)):
        ax.axhline(hval, color="0.35", lw=0.6, ls=(0, (1, 1.5)), zorder=1)
        ax.text(xmax, hval + 0.025*ytop, hlab, fontsize=LABEL_FS, color="0.25", ha="right", va="bottom")
        ax.plot(r1, sphere(chi)[k], color=C_PE, lw=1.4, zorder=3)
        ax.plot(chi, sphere(chi)[k], color=C_ELL, lw=1.4, zorder=3)
        ax.plot(h, sphere(chi_h)[k], color=C_HYP, lw=1.4, zorder=3)
        ax.plot(r, euclid[k], color=C_EUC, lw=1.0, ls=DASH, zorder=4)
        ax.plot([1.0], [sphere(np.pi/2)[k]], "o", ms=3.2, color=C_PE, mec=C_PE, mew=0.5, zorder=5)
        ax.plot([np.pi/2], [sphere(np.pi/2)[k]], "o", ms=3.2, color=C_ELL, mec=C_ELL, mew=0.5, zorder=5)
        ax.set_xlim(0, xmax); ax.set_ylim(0, ytop)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(ystep))
        ax.set_ylabel(ylab)
    axes[1].set_xlabel("radial coordinate")
    handles = [Line2D([], [], color=C_PE, lw=1.4), Line2D([], [], color=C_ELL, lw=1.4),
               Line2D([], [], color=C_HYP, lw=1.4), Line2D([], [], color=C_EUC, lw=1.0, ls=DASH)]
    labels = [r"pseudo-Euclidean, $r_1/\mathrm{R}$", r"elliptical, $\chi$",
              r"hyperbolic, $h/\mathrm{R}$", r"Euclidean, $r/\mathrm{R}$"]
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.53, -0.02), ncol=4,
               fontsize=LABEL_FS, handlelength=2.2, columnspacing=1.4)
    return fig


if __name__ == "__main__":
    print(f"mode {Z_MODE:.4f}, median {Z_MEDIAN:.4f}, fraction beyond z = 5: "
          f"{1 - n_below(5.0)/N_TOTAL:.3f}, total volume / R^3 {sphere(np.pi/2)[0]:.4f} (pi^2 = {np.pi**2:.4f})")
    save(cutoff_figure(), "cutoff_number_redshift")
    save(volume_figure(), "volume_area_circumference")
