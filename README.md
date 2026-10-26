# Magnitude–redshift data selection: CfA, SDSS and DESI DR1

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23587580.svg)](https://doi.org/10.5281/zenodo.23587580)

This repository reproduces the data selection, absolute-magnitude statistics and
comparison figures in T. Miller, *Radius, redshift, rotation, and radiation: the physical anti-de Sitter metric and entropic gravity* (2026), prepared for submission to JHEP. It compares the Hubble relation
(q₀ = 1) with the anti-de Sitter magnitude–redshift relation (appendix C of the paper) in three samples observed about a decade apart: the CfA Redshift Catalogue with the Hewitt & Burbidge quasar catalog (observations to 1995), SDSS-I/II (2000–2008) and DESI DR1 (2021–2022). The anti-de Sitter relation used here is eq. (7) of Miller & Miller (1995), Ap&SS 227, 53, which is called the de Sitter relation there.

| Script | Reproduces | Run time |
|---|---|---|
| `cfa_hb_selection.py` | CfA + Hewitt & Burbidge sample and its statistics, with the values of Miller & Miller (1995) for comparison | <1 min (includes a 7 MB download) |
| `sdss_legacy_selection.py` | SDSS-I/II sample and its statistics | ~10 min (includes a 1.7 GB download) |
| `desi_bright_selection.py` | DESI DR1 bright-time sample and its statistics | ~10 min (`--test`: 1% slice) |
| `make_figure.py` | the main figure (absolute magnitude; rows CfA, SDSS and DESI, columns anti-de Sitter and Hubble) and the appendix figure (apparent magnitude; CfA, SDSS and DESI); PDF and PNG | ~1 min |

## Quick start

```bash
pip install -r requirements.txt
python cfa_hb_selection.py             # writes cfa_hb_sample.csv.gz
python sdss_legacy_selection.py        # writes sdss_legacy_sample.csv.gz
python desi_bright_selection.py        # writes desi_dr1_main_bright_sample.csv.gz
python make_figure.py                  # writes absolute_magnitude_cfa_sdss_desi and apparent_magnitude_cfa_sdss_desi (.png/.pdf)
```

All data are public, and no accounts are needed. Each script prints the sample size and the statistics, so they can be checked against the tables below. Tested with Python 3.11.17, numpy 2.4.6, pandas 2.3.3, matplotlib 3.11.2, astropy 8.0.2, scipy 1.17.1 and requests 2.34.2.

## 1. CfA Redshift Catalogue and Hewitt & Burbidge quasars (to 1995)

This sample rebuilds the one in Miller & Miller (1995), Ap&SS 227, 53, from the latest versions of its catalogs, which are on the CDS archive (`https://cdsarc.cds.unistra.fr/ftp/cats/`):

- **Galaxies.** The CfA Redshift Catalogue, Version June 1995 (Huchra, Geller, Clemens, Tokarz & Michel 1995; CDS VII/193). The files are `zcat.dat` (57,536 entries, heliocentric velocities cz) and `zbig.dat` (1,202 high-redshift objects).
- **Quasars.** A Revised and Updated Catalog of Quasi-Stellar Objects (Hewitt & Burbidge 1993, ApJS 87, 451; CDS VII/158), file `table1_1.dat.gz` (7,315 records, complete to 1992 December 31).

The 1995 paper used the 1990 CfA catalog, ZBIG (Huchra & Clemens 1991) and Hewitt & Burbidge (1989), which are not available in machine-readable form.

| Step | Selection |
|---|---|
| galaxies, `zcat.dat` | cz > 0 and a magnitude, z = cz/c ≤ 0.33 (c = 299,792.5 km/s, as in the catalog); 26,754 |
| galaxies, `zbig.dat` | z > 0.33 and a magnitude; 527 |
| | CfA type −9 (QSO/AGN) excluded from both, because the quasars come from Hewitt & Burbidge |
| quasars | emission redshift z > 0 and a magnitude V; 7,031 (14 BL Lac objects) |
| duplicates | galaxies within 30″ and \|Δz\| < 0.02 of a Hewitt & Burbidge quasar are removed (73 zcat, 56 zbig) |
| quality cut | blue-band galaxy magnitudes only (below); removes 3,054 galaxies |

The sample contains **31,129** objects (24,098 galaxies and 7,031 quasars). The magnitude m is the CfA blue magnitude for galaxies and the Hewitt & Burbidge V magnitude (often a photographic magnitude) for quasars.

**Quality cut: blue-band magnitudes.** The CfA magnitude column is a blue magnitude, but some entries come from V, R or I photometry (source codes D, d, G, I, r, R, v, V, W) or have an unknown band (code 4, and blank codes in `zbig.dat`; in `zcat.dat` a blank code means B(0) of the First Reference Catalogue). Galaxies with these magnitudes are removed.

**Comparison with 1995.** The script prints the 1995 values next to its own. For all objects, the Hubble slope is −2.99 mag/dex with R² = 0.79 (1995: −3.0, R² = 0.81), and the anti-de Sitter slope is +0.14 mag/dex with R² = 0.008 (1995: +0.076, R² = 0.0026). Before the blue-band cut, these are −2.92 and +0.20.

**Flux-limited core.** The 14,750 galaxies with Zwicky magnitudes (source codes 1 and H), the flux-limited CfA survey proper, have median z = 0.021. On their own they give an anti-de Sitter slope of −1.76 mag/dex (R² = 0.29) and a Hubble slope of −4.32 (R² = 0.71). Within a flux-limited sample of this narrow redshift range, intrinsically brighter objects are seen to larger distances (Malmquist bias), which steepens both slopes. The script prints these values for reference.

## 2. SDSS-I/II (2000–2008)

**Spectra.** SDSS DR17 SkyServer (Abdurro'uf et al. 2022), `SpecObj` joined to `PhotoObj`, with `sciencePrimary = 1`, `zWarning = 0`, `z > 0`, class `GALAXY` or `QSO`, and `survey = 'sdss'`. This returns 936,278 spectra.

**Quasar redshifts.** These come from the DR16Q catalog (Lyke et al. 2020, ApJS 250, 8), file `DR16Q_v4.fits`, matched on (PLATE, MJD, FIBERID). The script checks the file's SHA-1 checksum after download.

- Matched spectra are quasars, with redshift `Z` from DR16Q.
- Unmatched `GALAXY` spectra are galaxies, with the SpecObj redshift.
- Unmatched `QSO` spectra are discarded.
- Finally, 0 < modelMag_r < 40 and redshift > 0 are required, and positive r-band flux, modelMag_r < 24.80 (quality cuts below).

The sample contains **921,067** objects (853,845 galaxies and 67,222 quasars).

## 3. DESI DR1, main survey, bright time (2021–2022)

The data are the DESI Data Release 1 redshift catalog (DESI Collaboration 2026, AJ 171, 285; arXiv:2503.14745), with the Legacy Surveys target photometry distributed with DR1. They are queried through the NOIRLab Astro Data Lab TAP service, `https://datalab.noirlab.edu/tap` (Nikutta et al. 2020).

```sql
SELECT z.targetid, z.z, z.spectype, p.flux_r, z.deltachi2
FROM desi_dr1.zpix AS z
JOIN desi_dr1.photometry AS p ON z.targetid = p.targetid
WHERE z.survey = 'main' AND z.program = 'bright'
  AND z.zcat_primary = 't' AND z.zwarn = 0
  AND z.spectype IN ('GALAXY','QSO')
  AND z.z > 0 AND p.flux_r > 0
```

| Cut | Purpose |
|---|---|
| `survey='main' AND program='bright'` | DESI main survey, bright-time program |
| `zcat_primary='t'` | one best spectrum per unique target |
| `zwarn=0` | DESI's flag for a redshift with no warnings |
| `spectype IN ('GALAXY','QSO')` | spectral class from the DESI pipeline (stars excluded) |
| `z>0`, `flux_r>0` | required to compute log z and m_r |

The query returns 6,015,420 objects. Galaxies must also pass the quality cuts below (`DELTACHI2 > 40` and z < 1.5); no other objects are removed. The sample contains **5,958,885** objects (5,890,529 galaxies and 68,356 quasars). The script splits the query into 100 equal slices of `random_id`, a uniform random column on [0, 100), purely to keep each download small.

## 4. Magnitudes and statistics

- m = m_r = 22.5 − 2.5 log₁₀(flux_r) for DESI (flux in nanomaggies, AB), m = modelMag_r for SDSS, and the catalog magnitudes for CfA (section 1). No extinction or K-correction is applied.
- Anti-de Sitter: M = m − 2.5 log₁₀[1 − (1+z)⁻²] + C_AdS, with C_AdS = −55.3
- Hubble (q₀ = 1): M = m − 5 log₁₀ z + C_H, with C_H = −43.0
- The constants were set for r-band magnitudes. A different band, or different constants, shifts every M by the same amount and does not change the widths or slopes.
- Δ68 = P84 − P16, Δ95 = P97.5 − P2.5, Δ99.7 = P99.85 − P0.15. These use `numpy.percentile` with linear interpolation, over all objects pooled.
- The slope is an unweighted least-squares fit of M on log₁₀ z, with no clipping.
- Running medians in the figures use bins of 0.1 in log₁₀ z and need at least 100 galaxies or 100 quasars per bin.
- Figure shading is the fraction of each sample per cell, on one logarithmic scale for all panels; cells holding a single object are a fixed light gray. Magnitudes increase upward in all panels.

| Sample | Relation | Δ68 | Δ95 | Δ99.7 [mag] | Slope [mag/dex] |
|---|---|---|---|---|---|
| CfA + Hewitt & Burbidge | anti-de Sitter | 2.3 | 5.2 | 8.8 | +0.1 |
| CfA + Hewitt & Burbidge | Hubble | 5.8 | 11.1 | 17.3 | −3.0 |
| CfA + Hewitt & Burbidge | apparent m (baseline) | 4.3 | 8.0 | 12.4 | +2.0 |
| SDSS-I/II | anti-de Sitter | 1.4 | 3.3 | 7.6 | +0.1 |
| SDSS-I/II | Hubble | 2.2 | 6.7 | 12.5 | −3.1 |
| SDSS-I/II | apparent m_r (baseline) | 2.1 | 4.7 | 8.2 | +1.9 |
| DESI DR1 bright | anti-de Sitter | 1.4 | 3.1 | 6.1 | −0.4 |
| DESI DR1 bright | Hubble | 2.0 | 4.5 | 10.8 | −3.4 |
| DESI DR1 bright | apparent m_r (baseline) | 1.6 | 3.5 | 5.7 | +1.6 |

All three samples are selected by apparent magnitude, so the spread of m provides a reference for the spread of M. The CfA sample combines catalogs with different magnitude limits and photometric systems; its spreads are wider than those of SDSS and DESI in all three rows.
For DESI, `desi_bright_selection.py` also prints the Hubble statistics for q₀ = 0.5 (2.0/4.6/11.2 mag)
and flat ΛCDM with Ωm = 0.3 (2.2/4.8/11.9 mag): smaller q₀ widens the Hubble distribution; of the cases considered,
q₀ = 1 gives the narrowest.

### Quality cuts

- **CfA: blue-band galaxy magnitudes** (section 1). This removes 3,054 galaxies.
- **SDSS: positive r-band flux.** SDSS magnitudes are asinh magnitudes (Lupton, Gunn & Szalay 1999), which stay finite at zero flux: in r, zero flux gives `modelMag_r` = 24.80 and negative flux gives fainter values. Requiring `modelMag_r` < 24.80 matches the DESI condition `flux_r > 0` and removes 504 objects (485 galaxies, 19 quasars), which would otherwise form faint diagonal sequences in the SDSS absolute-magnitude panels.
- **DESI galaxies: redshift quality.** Galaxies must satisfy DESI's redshift-success criterion for bright-time galaxies, `ZWARN = 0` and `DELTACHI2 > 40` (DESI Collaboration 2025, JCAP 07, 017), and z < 1.5 (Ratajczak et al. 2026, AJ 171, 71), which removes a pile-up below the z = 1.7 limit of the Redrock galaxy templates. Quasars are not cut. This removes 56,535 galaxies (0.94%): 51,735 with `DELTACHI2` ≤ 40 and 4,800 more at z ≥ 1.5.

All three selection scripts print the statistics before and after the cuts; the samples they write, and the figures, use the samples after the cuts.

| Before quality cuts | Relation | Δ68 | Δ95 | Δ99.7 [mag] | Slope [mag/dex] |
|---|---|---|---|---|---|
| CfA + Hewitt & Burbidge (34,183) | anti-de Sitter | 2.3 | 5.4 | 8.9 | +0.2 |
| CfA + Hewitt & Burbidge | Hubble | 5.3 | 10.9 | 16.9 | −2.9 |
| SDSS-I/II (921,571) | anti-de Sitter | 1.4 | 3.3 | 8.5 | +0.1 |
| SDSS-I/II | Hubble | 2.3 | 6.8 | 13.0 | −3.1 |
| DESI DR1 bright (6,015,420) | anti-de Sitter | 1.4 | 3.2 | 6.3 | −0.5 |
| DESI DR1 bright | Hubble | 2.0 | 4.5 | 11.5 | −3.5 |

## Citing this code

Please cite the Zenodo archive: T. Miller, *Magnitude–redshift analysis of SDSS and DESI galaxies and quasars*, Zenodo (2026), doi:[10.5281/zenodo.23587580](https://doi.org/10.5281/zenodo.23587580). This DOI always resolves to the latest version. Version v1.1.0 is doi:[10.5281/zenodo.23617979](https://doi.org/10.5281/zenodo.23617979); v1.0.0 is doi:[10.5281/zenodo.23587581](https://doi.org/10.5281/zenodo.23587581).

## Data credits

When using these data, please cite the CfA Redshift Catalogue (Huchra et al. 1995, CDS VII/193), the Hewitt & Burbidge (1993) quasar catalog, SDSS DR17 (Abdurro'uf et al. 2022), DR16Q (Lyke et al. 2020), DESI DR1 (DESI Collaboration 2026) and Astro Data Lab (Nikutta et al. 2020).

## License

MIT; see `LICENSE`.
