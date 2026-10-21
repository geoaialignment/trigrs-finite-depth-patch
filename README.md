# TRIGRS finite-depth solver: corrected kernel and pre-run exposure check

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23534840.svg)](https://doi.org/10.5281/zenodo.23534840) [![GeoAI Alignment](https://img.shields.io/badge/GeoAI%20Alignment-geoaialign.ai-1f6feb)](https://geoaialign.ai)

Maintained by **[GeoAI Alignment, Inc.](https://geoaialign.ai)**, developer of the G.AI.A GeoAI agent platform for spatial decision support.

This repository distributes a patch to [TRIGRS v2.1.0](https://github.com/usgs/landslides-trigrs/tree/v2.1.0) (Baum and Alvioli, 2016) that corrects the finite-depth pressure-head solver, together with a regression test against an independent analytic solution and a screen that tells a user, before a run, whether an input set is exposed to the departures.

In TRIGRS v2.1.0, results of the finite-depth solver depart from the analytic solution of the same equation once a forcing period becomes older than a dimensionless age A = D(t − a)/L² of 5, where D is the slope-corrected diffusivity, L the soil depth and t − a the time since the period began. No tested case below that age had a material error, so many short-event applications are unaffected. Longer simulations, such as runs with an antecedent dry or wet lead, can be exposed. Once older periods reach an age of about 72 to 76, later rainfall can be dropped altogether.

## What is here

| Path | Contents |
|---|---|
| `release_package/release.patch` | Patch to ten Fortran files under `src/TRIGRS/` of TRIGRS v2.1.0 |
| `release_package/build.py` | Build of the serial and MPI executables with gfortran and MPICH |
| `release_package/regression/run_regression.py` | Six-case regression against an independent analytic reference |
| `release_package/screen/screen.py` | Pre-run exposure check of a `tr_in.txt` input set; no solver run needed |
| `release_package/logs/` | Build, patch, regression and screen logs from the verification runs |

Full instructions are in [`release_package/README.md`](release_package/README.md).

## Quick start

```sh
git clone https://github.com/geoaialignment/trigrs-finite-depth-patch.git
cd trigrs-finite-depth-patch
git clone --branch v2.1.0 --single-branch https://github.com/usgs/landslides-trigrs.git
(cd landslides-trigrs/src/TRIGRS && patch -p0 --batch --fuzz=0 -i ../../../release_package/release.patch)
python3 release_package/build.py landslides-trigrs/src/TRIGRS build_serial
python3 release_package/regression/run_regression.py build_serial/trg --workdir regression_run
```

To check whether your own input set is exposed before running TRIGRS:

```sh
python3 release_package/screen/screen.py path/to/tr_in.txt
```

## Citation

The analysis is described in a manuscript submitted to Geoscientific Model Development:

> Park, E., Kim, T., and Park, J.: Analytic verification of the finite-depth pressure-head solver of TRIGRS (versions 2.0.06b and 2.1.0), with a pre-run exposure check and a corrected kernel, submitted to Geosci. Model Dev., 2026.

Code, data and the corrected kernel are archived at Zenodo:

> Park, E., Kim, T., and Park, J.: Code, data and corrected kernel for "Analytic verification of the finite-depth pressure-head solver of TRIGRS (versions 2.0.06b and 2.1.0), with a pre-run exposure check and a corrected kernel" (v1.0.0), Zenodo, https://doi.org/10.5281/zenodo.23534840, 2026.

## Licence

The patch, scripts and documentation are dedicated to the public domain under [CC0 1.0](LICENSE.md), consistent with the public-domain status of the upstream USGS code. This is not an official USGS release. We thank the TRIGRS developers, Rex L. Baum and Massimiliano Alvioli, for making the source and its revision notes public.

## GeoAI Alignment

[GeoAI Alignment, Inc.](https://geoaialign.ai) develops **G.AI.A**, a GeoAI agent platform that connects spatial data, satellite and field observations with expert analysis and simulation tools. A user states the area and the question in plain language; the agents run the analyses and return the results with the evidence behind them.

Applications include landslide, flood and wildfire risk, groundwater and surface water, agriculture and forestry, ecology, urban planning and energy siting. See [solutions](https://geoaialign.ai/solutions) and our [research papers](https://geoaialign.ai/blog?category=research#papers).

For landslide hazard analysis, model verification or a G.AI.A demonstration for your region, [contact us](https://geoaialign.ai/contact) or write to contact@geoaialign.com.

This work was carried out by GeoAI Alignment, Inc. with the School of Earth System Sciences, Kyungpook National University, Daegu, Republic of Korea. Scientific contact: Eungyu Park (egpark@knu.ac.kr).
