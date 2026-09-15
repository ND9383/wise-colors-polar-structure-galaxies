# Mid-Infrared Colours of Polar-Structure Galaxies in COUGS-DESI

Analysis code for:

**Dufera, N. (2026). "Mid-Infrared Colours Trace Formation Channel in the
COUGS-DESI Catalogue of Polar-Structure Galaxies." MNRAS Letters
(submitted).**

## Overview

This repository contains the full analysis pipeline used to cross-match
the 2,983 usable polar-structure galaxies (PSGs) in the COUGS-DESI
catalogue (Bahr-Kalus & Mosenkov 2026) against the AllWISE mid-infrared
catalogue, and to test whether the catalogue's structural subtypes differ
in WISE colour in the direction predicted by their proposed formation
channels.

## Data sources

Both datasets are public and are pulled automatically via `astroquery`:

- **COUGS-DESI**: VizieR catalogue `J/A+A/710/A145`
- **AllWISE**: VizieR catalogue `II/328/allwise` (Cutri et al. 2013)

No manual downloads or authentication are required.

## Requirements

```
pip install astroquery astropy scipy pandas numpy
```

## Usage

```
python wise_analysis.py
```

This reproduces every statistical result reported in the Letter, in the
order they appear in the paper: the omnibus subtype comparison
(Kruskal-Wallis), the redshift-dependence check, the PR-vs-PTS core
prediction test, the PB quiescent-fraction and star-forming-fraction
comparisons using the literature-correct three-zone WISE classification
(Jarrett et al. 2011, 2017), and the AGN fraction test using the Stern
et al. (2012) criterion.

**Note:** the full cross-match queries AllWISE once per galaxy (2,983
individual queries) and takes approximately 30-40 minutes to run given
typical VizieR response times. This is expected, not an error.

## Method summary

WISE colours (W1$-$W2, W2$-$W3) are computed from the nearest AllWISE
counterpart within 6 arcsec of each COUGS-DESI target. Galaxies are
classified using the standard three-zone scheme: quiescent (W2$-$W3 <
2.0), intermediate (2.0-3.5), and star-forming (> 3.5), with AGN
identified separately via W1$-$W2 > 0.8 restricted to W2 < 15.05 mag
(Stern et al. 2012).

## Citation

If you use this code, please cite the associated Letter (citation details
to be added upon publication).
