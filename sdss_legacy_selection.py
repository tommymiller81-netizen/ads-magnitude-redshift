"""
Reproduce the SDSS legacy (SDSS-I/II) comparison sample and its statistics.

Spectra and photometry: SDSS DR17 SkyServer (Abdurro'uf et al. 2022).
Quasar redshifts:       DR16Q quasar catalog, DR16Q_v4.fits (Lyke et al. 2020),
                        downloaded from the SDSS Science Archive Server (1.7 GB).

Selection:
  SpecObj: sciencePrimary = 1, zWarning = 0, z > 0, class IN ('GALAXY','QSO'), survey = 'sdss'
  Join DR16Q on (PLATE, MJD, FIBERID):
    matched              -> quasar, redshift = DR16Q Z
    unmatched, GALAXY    -> galaxy, redshift = SpecObj z
    unmatched, QSO       -> discarded
  Then 0 < modelMag_r < 40 and redshift > 0 (redshifts kept at full precision).
Quantities and statistics as in desi_bright_selection.py, with m_r = modelMag_r.

Optional zero-flux test (--exclude-zero-flux; appendix C of the paper):
  SDSS magnitudes are asinh magnitudes (Lupton, Gunn & Szalay 1999). In r, zero flux gives
  modelMag_r = 24.80 (softening b_r = 1.2e-10) and negative flux gives fainter values. With the
  flag, the statistics are also printed without the objects at modelMag_r >= 24.80. The saved
  sample always contains all objects.

Usage:  python sdss_legacy_selection.py
        python sdss_legacy_selection.py --exclude-zero-flux   # also run the zero-flux test
Requires: numpy, pandas, requests, astropy.
"""
import io, os, sys, time, hashlib, requests, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits

SKY = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
DR16Q_URL = "https://data.sdss.org/sas/dr16/eboss/qso/DR16Q/DR16Q_v4.fits"
DR16Q_SHA1 = "9c9db93f1fb23b31576dd04836ede78ff982f9cb"
ZERO_FLUX = "--exclude-zero-flux" in sys.argv
M0_R = -2.5/np.log(10)*np.log(1.2e-10)     # asinh r magnitude at zero flux, 24.802
SQL = ("SELECT s.plate, s.mjd, s.fiberID, s.z, s.class, p.modelMag_r "
       "FROM SpecObj AS s JOIN PhotoObj AS p ON s.bestObjID = p.objID "
       "WHERE s.sciencePrimary = 1 AND s.zWarning = 0 AND s.z > 0 "
       "AND s.class IN ('GALAXY','QSO') AND s.survey = 'sdss' AND (s.fiberID % 40) = {k}")

def chunk(k):
    # the fiberID modulo only splits the query into 40 parts to stay under SkyServer's row limit
    for attempt in range(4):
        t = requests.get(SKY, params={"cmd": SQL.format(k=k), "format": "csv"}, timeout=900).text
        if t.startswith("#Table"): t = t.split("\n", 1)[1]
        if t.startswith("plate"): return pd.read_csv(io.StringIO(t))
        time.sleep(10)
    raise RuntimeError(t[:300])

def get_dr16q(fn="DR16Q_v4.fits"):
    if not os.path.exists(fn):
        with requests.get(DR16Q_URL, stream=True, timeout=3600) as r, open(fn + ".part", "wb") as f:
            r.raise_for_status()
            for b in r.iter_content(1 << 22): f.write(b)
        os.replace(fn + ".part", fn)
    h = hashlib.sha1()
    with open(fn, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    assert h.hexdigest() == DR16Q_SHA1, "DR16Q_v4.fits checksum mismatch"
    with fits.open(fn, memmap=True) as hd:
        d = hd[1].data
        return pd.DataFrame({"plate": np.asarray(d["PLATE"]).astype(int), "mjd": np.asarray(d["MJD"]).astype(int),
                             "fiberID": np.asarray(d["FIBERID"]).astype(int), "z16": np.asarray(d["Z"]).astype(float)})

def main():
    with ThreadPoolExecutor(4) as ex:
        s = pd.concat(list(ex.map(chunk, range(40))), ignore_index=True)
    print(f"SpecObj query: {len(s):,} spectra")
    s = s.merge(get_dr16q(), on=["plate", "mjd", "fiberID"], how="left")
    s["sample"] = np.where(s.z16.notna(), "QSO", np.where(s["class"] == "GALAXY", "GALAXY", "drop"))
    s["z_use"] = np.where(s.z16.notna(), s.z16, s.z)
    k = s[(s["sample"] != "drop") & (s.z_use > 0) & (s.modelMag_r > 0) & (s.modelMag_r < 40)]
    k.to_csv("sdss_legacy_sample.csv.gz", index=False)
    report(k)
    if ZERO_FLUX:
        zf = k.modelMag_r >= M0_R - 5e-4
        print(f"\nZero-flux test: removing {zf.sum():,} objects with modelMag_r >= {M0_R:.2f} "
              f"({(k[zf]['sample']=='GALAXY').sum():,} galaxies, {(k[zf]['sample']=='QSO').sum():,} quasars)")
        report(k[~zf])

def report(k):
    lz = np.log10(k.z_use.to_numpy()); m = k.modelMag_r.to_numpy()
    M = {"anti-de Sitter (C_AdS = -55.3)": m - 2.5*np.log10(1 - (1 + k.z_use.to_numpy())**-2) - 55.3,
         "Hubble (q0 = 1, C_H = -43.0)":   m - 5*lz - 43.0,
         "apparent m_r (no distance)":     m}
    print(f"N = {len(k):,}  ({(k['sample']=='GALAXY').sum():,} galaxies, {(k['sample']=='QSO').sum():,} quasars)")
    for name, v in M.items():
        p = np.percentile(v, [0.15, 2.5, 16, 84, 97.5, 99.85]); slope = np.polyfit(lz, v, 1)[0]
        print(f"{name}: D68 = {p[3]-p[2]:.1f}  D95 = {p[4]-p[1]:.1f}  D99.7 = {p[5]-p[0]:.1f} mag;"
              f"  LSQ slope = {slope:.1f} mag/dex")

if __name__ == "__main__":
    main()
