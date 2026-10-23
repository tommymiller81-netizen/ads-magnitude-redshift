# Magnitude–redshift data selection: DESI DR1 and SDSS

This repository reproduces the data selection, absolute-magnitude statistics and
comparison figure in T. Miller, *Radius, redshift, rotation, and radiation: the physical anti-de Sitter metric and entropic gravity* (2026), prepared for submission to JHEP. It compares the Hubble relation
(q₀ = 1) with the anti-de Sitter magnitude–redshift relation (eq. 21 of the paper).

| Script | Reproduces | Run time |
|---|---|---|
| `desi_bright_selection.py` | DESI DR1 bright-time sample and its statistics | ~10 min (`--test`: 1% slice) |
| `make_figure.py` | the comparison figure (PNG and PDF) | ~1 min |
| `sdss_legacy_selection.py` | SDSS legacy comparison sample and its statistics | ~10 min (includes a 1.7 GB download) |

## Quick start

```bash
pip install -r requirements.txt
python desi_bright_selection.py        # writes desi_dr1_main_bright_sample.csv.gz
python make_figure.py                  # writes desi_bright_hubble_vs_ads.png/.pdf
python sdss_legacy_selection.py        # writes sdss_legacy_sample.csv.gz
```

All data are public, and no accounts are needed. Each script prints the sample size and the statistics, so they can be checked against the tables below. Tested with Python 3.11.17, numpy 2.4.6, pandas 2.3.3, matplotlib 3.11.2, astropy 8.0.2 and requests 2.34.2.

## 1. DESI DR1, main survey, bright time

The data are the DESI Data Release 1 redshift catalogue (DESI Collaboration 2026, AJ 171, 285; arXiv:2503.14745), with the Legacy Surveys target photometry distributed with DR1. They are queried through the NOIRLab Astro Data Lab TAP service, `https://datalab.noirlab.edu/tap` (Nikutta et al. 2020).

```sql
SELECT z.targetid, z.z, z.spectype, p.flux_r
FROM desi_dr1.zpix AS z
JOIN desi_dr1.photometry AS p ON z.targetid = p.targetid
WHERE z.survey = 'main' AND z.program = 'bright'
  AND z.zcat_primary = 't' AND z.zwarn = 0
  AND z.spectype IN ('GALAXY','QSO')
  AND z.z > 0 AND p.flux_r > 0
```

| Cut | Purpose |
|---|---|
| `survey='main' AND program='bright'` | DESI main survey, bright-time programme |
| `zcat_primary='t'` | one best spectrum per unique target |
| `zwarn=0` | DESI's flag for a redshift with no warnings |
| `spectype IN ('GALAXY','QSO')` | spectral class from the DESI pipeline (stars excluded) |
| `z>0`, `flux_r>0` | required to compute log z and m_r |

No other objects are removed; known redshift outliers are kept. The sample contains **6,015,420** objects (5,947,064 galaxies and 68,356 quasars). The script splits the query into 100 equal slices of `random_id`, a uniform random column on [0, 100), purely to keep each download small.

## 2. SDSS legacy (comparison sample)

**Spectra.** SDSS DR17 SkyServer (Abdurro'uf et al. 2022), `SpecObj` joined to `PhotoObj`, with `sciencePrimary = 1`, `zWarning = 0`, `z > 0`, class `GALAXY` or `QSO`, and `survey = 'sdss'`. This returns 936,278 spectra.

**Quasar redshifts.** These come from the DR16Q catalogue (Lyke et al. 2020, ApJS 250, 8), file `DR16Q_v4.fits`, matched on (PLATE, MJD, FIBERID). The script checks the file's SHA-1 checksum after download.

- Matched spectra are quasars, with redshift `Z` from DR16Q.
- Unmatched `GALAXY` spectra are galaxies, with the SpecObj redshift.
- Unmatched `QSO` spectra are discarded.
- Finally, 0 < modelMag_r < 40 and redshift > 0 are required.

The sample contains **921,571** objects (854,330 galaxies and 67,241 quasars).

## 3. Magnitudes and statistics

- m_r = 22.5 − 2.5 log₁₀(flux_r) for DESI (flux in nanomaggies, AB), and m_r = modelMag_r for SDSS. No extinction or K-correction is applied.
- Hubble (q₀ = 1): M = m_r − 5 log₁₀ z − 43.0
- Anti-de Sitter (eq. 21): M = m_r − 2.5 log₁₀[1 − (1+z)⁻²] − 55.3
- Δ68 = P84 − P16, Δ95 = P97.5 − P2.5, Δ99.7 = P99.85 − P0.15. These use `numpy.percentile` with linear interpolation, over all objects pooled.
- The slope is an unweighted least-squares fit of M on log₁₀ z, with no clipping.
- Running medians in the figure use bins of 0.1 in log₁₀ z, plotted only when a bin holds at least 1000 galaxies, or at least 200 quasars.

| Sample | Relation | Δ68 | Δ95 | Δ99.7 [mag] | Slope [mag/dex] |
|---|---|---|---|---|---|
| DESI DR1 bright | Hubble | 1.97 | 4.53 | 11.5 | −3.46 |
| DESI DR1 bright | anti-de Sitter | 1.39 | 3.15 | 6.30 | −0.471 |
| SDSS legacy | Hubble | 2.25 | 6.76 | 13.0 | −3.10 |
| SDSS legacy | anti-de Sitter | 1.38 | 3.30 | 8.51 | +0.0538 |

## Data credits

When using these data, please cite DESI DR1 (DESI Collaboration 2026), SDSS DR17 (Abdurro'uf et al. 2022), DR16Q (Lyke et al. 2020) and Astro Data Lab (Nikutta et al. 2020).

## Licence

MIT; see `LICENSE`.
