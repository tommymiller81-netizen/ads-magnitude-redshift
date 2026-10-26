"""
Reproduce the CfA + Hewitt & Burbidge comparison sample (the sample of Miller & Miller 1995,
Ap&SS 227, 53, rebuilt from the latest versions of its catalogs) and its statistics.

Data (CDS archive, no account needed):
  CfA Redshift Catalogue, Version June 1995 (Huchra, Geller, Clemens, Tokarz & Michel 1995;
    CDS VII/193): zcat.dat (57,536 entries, heliocentric velocities cz in km/s) and zbig.dat
    (1,202 high-redshift objects, redshifts z).
  A Revised and Updated Catalog of Quasi-Stellar Objects (Hewitt & Burbidge 1993, ApJS 87, 451;
    CDS VII/158): table1_1.dat.gz (7,315 records), complete to 1992 December 31.
The 1995 paper used the 1990 CfA catalog, ZBIG (1991) and Hewitt & Burbidge (1989), which are
not available in machine-readable form; the versions above are their direct successors.

Selection, as in the 1995 paper:
  galaxies: zcat entries with Vh > 0 and a magnitude, z = Vh/c <= 0.33 (c = 299792.5 km/s, as in
            the catalog), plus zbig entries with z > 0.33 and a magnitude; CfA type -9 (QSO/AGN)
            excluded, because quasars come from Hewitt & Burbidge
  quasars:  Hewitt & Burbidge objects with an emission redshift z > 0 and a magnitude V
  duplicates: a galaxy within 30 arcsec and |dz| < 0.02 of a Hewitt & Burbidge quasar is removed
Quality cut: blue-band galaxy magnitudes only. The CfA magnitude column is a B magnitude, but some
  entries come from V, R or I photometry (source codes D, d, G, I, r, R, v, V, W) or have an
  unknown band (code 4, and blank codes in zbig; in zcat a blank code means B(0) of the First
  Reference Catalogue). Quasars keep the Hewitt & Burbidge V magnitude (often only m_pg).
  Statistics are printed before and after the cut; the sample written to cfa_hb_sample.csv.gz
  (and used for the figures) is the one after the cut.
Quantities and statistics as in desi_bright_selection.py, with m = the catalog magnitude, plus the
  squared correlation coefficient R^2 reported in the 1995 paper. The constants C_AdS and C_H were
  set for r-band magnitudes; a different band only shifts every M by a constant and does not
  change the widths or slopes.

Usage:  python cfa_hb_selection.py
Requires: numpy, pandas, requests, scipy.
"""
import gzip, os, requests, numpy as np, pandas as pd
from scipy.spatial import cKDTree

CDS = "https://cdsarc.cds.unistra.fr/ftp/cats/"
C = 299792.5
ZMAX_ZCAT = 0.33
BLUE = set("0123456789ABCHJLP") - {"4"}            # blue-band magnitude source codes (ReadMe of VII/193)
ZWICKY = {"1", "H"}                                 # Zwicky-catalog magnitudes (the CfA survey system)

ZCAT = [("name", 0, 11), ("RAh", 11, 13), ("RAm", 13, 15), ("RAs", 15, 19), ("DEsign", 19, 20),
        ("DEd", 20, 22), ("DEm", 22, 24), ("DEs", 24, 26), ("mag", 26, 31), ("Vh", 31, 36),
        ("mag_src", 39, 40), ("zRef", 40, 44), ("Ttype", 44, 46)]
ZBIG = [("name", 0, 11), ("RAh", 11, 13), ("RAm", 13, 15), ("RAs", 15, 19), ("DEsign", 19, 20),
        ("DEd", 20, 22), ("DEm", 22, 24), ("DEs", 24, 26), ("z", 27, 33), ("mag", 34, 39),
        ("mag_src", 40, 41), ("zRef", 42, 46), ("Ttype", 46, 48)]
HB = [("name", 0, 8), ("RAh", 9, 11), ("RAm", 12, 14), ("RAs", 15, 20), ("DEsign", 21, 22),
      ("DEd", 22, 24), ("DEm", 25, 27), ("DEs", 28, 32), ("BLLac", 33, 34), ("mag", 76, 81), ("z", 95, 100)]
TEXT = ("name", "DEsign", "mag_src", "BLLac")

# Miller & Miller (1995), Table I and figure 2: (slope, R^2)
REF = {"galaxies": {"Hubble": (-3.0, 0.60), "anti-de Sitter": (-0.27, 0.011)},
       "quasars":  {"Hubble": (-3.6, 0.52), "anti-de Sitter": (0.33, 0.0089)},
       "all":      {"Hubble": (-3.0, 0.81), "anti-de Sitter": (0.076, 0.0026), "apparent m": (2.0, 0.65)}}

def get(path):
    fn = os.path.basename(path)
    if not os.path.exists(fn):
        r = requests.get(CDS + path, timeout=300); r.raise_for_status()
        open(fn, "wb").write(r.content)
    return fn

def read(path, cols):
    fn = get(path)
    src = gzip.open(fn, "rt") if fn.endswith(".gz") else fn
    d = pd.read_fwf(src, colspecs=[(a, b) for _, a, b in cols], names=[c for c, _, _ in cols],
                    dtype=str, keep_default_na=False)
    for c in d.columns:
        d[c] = d[c].str.strip() if c in TEXT else pd.to_numeric(d[c].str.strip().replace("", np.nan), errors="coerce")
    d["ra1950"] = 15*(d.RAh + d.RAm/60 + d.RAs.fillna(0)/3600)
    d["dec1950"] = np.where(d.DEsign == "-", -1, 1)*(d.DEd + d.DEm/60 + d.DEs.fillna(0)/3600)
    return d

def unit(ra, dec):
    r, d = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.c_[np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)]

def report(k, title, ref=True):
    print(f"\n{title}")
    groups = [("galaxies", k[k["sample"] == "GALAXY"]), ("quasars", k[k["sample"] == "QSO"]), ("all", k)]
    if k["sample"].nunique() == 1:
        groups = groups[2:]
    for g, d in groups:
        if len(d) == 0:
            continue
        z = d.z.to_numpy(float); lz = np.log10(z); m = d.mag.to_numpy(float)
        M = {"anti-de Sitter": m - 2.5*np.log10(1 - (1 + z)**-2) - 55.3,
             "Hubble": m - 5*lz - 43.0, "apparent m": m}
        print(f"  {g}: N = {len(d):,}")
        for name, v in M.items():
            p = np.percentile(v, [0.15, 2.5, 16, 84, 97.5, 99.85])
            slope = np.polyfit(lz, v, 1)[0]; r2 = np.corrcoef(lz, v)[0, 1]**2
            r95 = REF.get(g, {}).get(name) if ref else None
            r95 = f"   (1995: slope {r95[0]:+}, R^2 {r95[1]})" if r95 else ""
            print(f"    {name:15s} D68 = {p[3]-p[2]:4.1f}  D95 = {p[4]-p[1]:4.1f}  D99.7 = {p[5]-p[0]:4.1f} mag;"
                  f"  slope = {slope:+.2f} mag/dex, R^2 = {r2:.4f}{r95}")

def main():
    zc = read("VII/193/zcat.dat", ZCAT); zb = read("VII/193/zbig.dat", ZBIG)
    hb = read("VII/158/table1_1.dat.gz", HB)
    print(f"zcat.dat {len(zc):,} entries; zbig.dat {len(zb):,}; Hewitt & Burbidge {len(hb):,} records")
    q = hb[(hb.z > 0) & hb.mag.notna()]
    zc = zc[(zc.Vh > 0) & zc.mag.notna() & (zc.Ttype != -9)].assign(z=lambda d: d.Vh/C, catalog="zcat")
    zc = zc[zc.z <= ZMAX_ZCAT]
    zb = zb[(zb.z > ZMAX_ZCAT) & zb.mag.notna() & (zb.Ttype != -9)].assign(catalog="zbig")
    print(f"galaxies: zcat (z <= {ZMAX_ZCAT}) {len(zc):,}, zbig (z > {ZMAX_ZCAT}) {len(zb):,};"
          f" quasars with z > 0 and V: {len(q):,} ({(q.BLLac == 'B').sum()} BL Lac objects)")
    tq = cKDTree(unit(q.ra1950, q.dec1950)); keep = []
    for name, g in (("zcat", zc), ("zbig", zb)):
        dist, idx = tq.query(unit(g.ra1950, g.dec1950), distance_upper_bound=np.radians(30/3600))
        hit = np.isfinite(dist); dup = np.zeros(len(g), bool)
        dup[hit] = np.abs(g.z.to_numpy()[hit] - q.z.to_numpy()[idx[hit]]) < 0.02
        print(f"{name} galaxies within 30 arcsec and |dz| < 0.02 of a Hewitt & Burbidge quasar (removed): {dup.sum()}")
        keep.append(g[~dup])
    cols = ["name", "ra1950", "dec1950", "z", "sample", "mag", "mag_src", "catalog"]
    k = pd.concat([*(g.assign(sample="GALAXY") for g in keep),
                   q.assign(sample="QSO", mag_src="HB", catalog="hewitt_burbidge")], ignore_index=True)[cols]
    k["ra1950"] = k.ra1950.round(5); k["dec1950"] = k.dec1950.round(5)
    report(k, "Before quality cut (all magnitude bands):")
    blue = (k["sample"] == "QSO") | k.mag_src.isin(BLUE) | ((k.mag_src == "") & (k.catalog == "zcat"))
    print(f"\nQuality cut removes {(~blue).sum():,} galaxies without a blue-band magnitude")
    k = k[blue]
    k.rename(columns={"mag": "Bmag"}).to_csv("cfa_hb_sample.csv.gz", index=False)
    report(k, f"Sample (after quality cut), written to cfa_hb_sample.csv.gz: {len(k):,} objects,"
              f" {(k['sample'] == 'QSO').mean():.0%} quasars")
    zw = k[(k["sample"] == "GALAXY") & k.mag_src.isin(ZWICKY)]
    report(zw, f"For reference, galaxies with Zwicky magnitudes only"
           f" (flux-limited CfA survey core; median z = {zw.z.median():.4f}):", ref=False)

if __name__ == "__main__":
    main()
