"""
Reproduce the DESI DR1 bright-time sample and the magnitude statistics.

Data: DESI Data Release 1 (DESI Collaboration 2025, arXiv:2503.14745), queried
through the NOIRLab Astro Data Lab TAP service (tables desi_dr1.zpix and
desi_dr1.photometry). No account is needed for these public tables.

Selection (all applied server-side):
  survey = 'main' AND program = 'bright'     DESI main survey, bright-time programme
  zcat_primary = 't'                         one best spectrum per TARGETID
  zwarn = 0                                  no redshift warning flags
  spectype IN ('GALAXY','QSO')               Redrock spectral class
  z > 0 AND flux_r > 0                       needed for log z and m_r
Nothing else is removed.

Quantities:
  m_r   = 22.5 - 2.5 log10(flux_r)          flux_r in nanomaggies (AB), no extinction
                                             or K-correction applied
  M_H   = m_r - 5 log10(z) - 43.0            Hubble, q0 = 1
  M_AdS = m_r - 2.5 log10[1 - (1+z)^-2] - 55.3
Band widths: D68 = P84 - P16, D95 = P97.5 - P2.5, D99.7 = P99.85 - P0.15
(numpy.percentile, default linear interpolation, all objects pooled).
Slope: ordinary least squares of M on log10 z (numpy.polyfit, degree 1).

Optional redshift-failure test (--exclude-failures; appendix C of the paper):
  also downloads DELTACHI2 (chi^2 difference between the best and second-best
  Redrock fits) and repeats the statistics after removing
    DELTACHI2 <= 40                            low-confidence redshifts
    spectype = 'GALAXY' AND z >= 1.6           pile-up below the z = 1.7 limit of
                                               the Redrock galaxy templates
  The full sample is still written and used for the figure; the reduced sample
  is written separately to desi_dr1_main_bright_sample_nofail.csv.gz.

Usage:  python desi_bright_selection.py                      # full sample (~10 min)
        python desi_bright_selection.py --test               # 1% slice, quick check
        python desi_bright_selection.py --exclude-failures   # also run the failure test
Requires: python >= 3.9, numpy, pandas, requests.
"""
import sys, io, time, requests, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor

TAP = "https://datalab.noirlab.edu/tap/sync"
WHERE = ("z.survey = 'main' AND z.program = 'bright' AND z.zcat_primary = 't' "
         "AND z.zwarn = 0 AND z.spectype IN ('GALAXY','QSO') AND z.z > 0 AND p.flux_r > 0")
COLS = "z.targetid, z.z, z.spectype, p.flux_r"
FAIL = "--exclude-failures" in sys.argv
if FAIL:
    COLS += ", z.deltachi2"

def query(adql):
    for attempt in range(4):
        r = requests.post(TAP, data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                     "MAXREC": "5000000", "QUERY": adql}, timeout=3600)
        if r.status_code == 200 and r.text.startswith("targetid"):
            return pd.read_csv(io.StringIO(r.text))
        time.sleep(20)
    raise RuntimeError(r.text[:300])

def dl_lcdm(z, om=0.3):
    """Flat LCDM luminosity distance in units of c/H0 (trapezoidal integral on a fine grid)."""
    zg = np.linspace(0, z.max()*1.001, 200001)
    dc = np.concatenate([[0], np.cumsum(np.diff(zg)*0.5*(1/np.sqrt(om*(1+zg[1:])**3 + 1-om)
                                                          + 1/np.sqrt(om*(1+zg[:-1])**3 + 1-om)))])
    return (1 + z)*np.interp(z, zg, dc)

def stats(d, label):
    m = 22.5 - 2.5*np.log10(d.flux_r.astype(float))
    z = d.z.astype(float); lz = np.log10(z)
    M = {"anti-de Sitter (C_AdS = -55.3)": m - 2.5*np.log10(1 - (1 + z)**-2) - 55.3,
         "Hubble (q0 = 1, C_H = -43.0)":   m - 5*lz - 43.0,
         # reference rows: smaller q0 and flat LCDM (constants omitted; widths do not depend on them)
         "Hubble (q0 = 0.5, Mattig)":      m - 5*np.log10(2*(z + 1 - np.sqrt(1 + z))),
         "Hubble (flat LCDM, Om = 0.3)":   m - 5*np.log10(dl_lcdm(z.to_numpy())),
         "apparent m_r (no distance)":     m}
    print(f"{label}: N = {len(d):,}  ({(d.spectype=='GALAXY').sum():,} galaxies, "
          f"{(d.spectype=='QSO').sum():,} quasars)")
    for name, v in M.items():
        p = np.percentile(v, [0.15, 2.5, 16, 84, 97.5, 99.85])
        slope = np.polyfit(lz, v, 1)[0]
        print(f"  {name}: D68 = {p[3]-p[2]:.1f}  D95 = {p[4]-p[1]:.1f}  D99.7 = {p[5]-p[0]:.1f} mag;"
              f"  LSQ slope = {slope:.1f} mag/dex")

def chunk(k):
    # random_id is uniform on [0, 100); slicing by it splits the query into 100 equal parts
    hi = "" if k == 99 else f" AND z.random_id < {k+1}"
    return query(f"SELECT {COLS} FROM desi_dr1.zpix AS z JOIN desi_dr1.photometry AS p "
                 f"ON z.targetid = p.targetid WHERE {WHERE} AND z.random_id >= {k}{hi}")

def main():
    ks = [0] if "--test" in sys.argv else range(100)
    with ThreadPoolExecutor(4) as ex:
        d = pd.concat(list(ex.map(chunk, ks)), ignore_index=True)
    d.to_csv("desi_dr1_main_bright_sample.csv.gz", index=False)
    stats(d, "Full sample")
    if FAIL:
        bad = (d.deltachi2 <= 40) | ((d.spectype == "GALAXY") & (d.z >= 1.6))
        print(f"Removed as likely redshift failures: {bad.sum():,} ({100*bad.mean():.2f}%)")
        d[~bad].to_csv("desi_dr1_main_bright_sample_nofail.csv.gz", index=False)
        stats(d[~bad], "Without likely failures")

if __name__ == "__main__":
    main()
