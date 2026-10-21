# TRIGRS v2.1.0 corrected kernel

Public software package accompanying the GMD study. The patch changes ten
Fortran routines and adds no instrumentation, input fields or output formats.
It preserves the existing pressure cap and factor-of-safety conventions.

## Corrections

**Complete finite-layer response.** In `savage`, `svgstp` and `svgstp_p`, the
late-time finite-layer response needs its mean (water-storage) and polynomial
terms as well as the decaying eigenmodes. The corrected step kernel retains
all these terms when convolving each rainfall period. This restores the
response to a uniform storm and redistribution after a burst, including
periods that enter the late-time branch.

**History-loop early exit.** A period whose two response tails evaluate to
zero does not establish that all subsequent periods contribute zero. This
can occur for an old dry lead while a younger storm still contributes. The
patch removes the early exit from the forcing-history loop in all ten
routines, allowing every relevant period to contribute.

**Transient-pressure assignment.** In `ivestp`, `ivestp_p`, `svgstp` and
`svgstp_p`, transient pressure must be assigned at each evaluated time even
when the newly computed response is zero. Removing the nonzero-response
assignment guard prevents an earlier value from remaining in the transient
pressure array. The finite-layer kernel also supplies the nonzero stored
water response needed during dry periods.

Changed files, all under upstream `src/TRIGRS/`:

- `iverson.f95`
- `ivestp.f95`
- `ivestp_p.f95`
- `pstpf.f95`
- `pstpf_p.f95`
- `pstpi.f95`
- `pstpi_p.f95`
- `savage.f95`
- `svgstp.f95`
- `svgstp_p.f95`

## Obtain and apply

Obtain the [USGS TRIGRS repository at tag v2.1.0](https://github.com/usgs/landslides-trigrs/tree/v2.1.0).
The pinned commit is `969c409f7bc63e7616ebe34dfd5f6afd967f2d2d`.
Run these commands from a directory containing this package as `release_package/`:

```sh
git clone --branch v2.1.0 --single-branch https://github.com/usgs/landslides-trigrs.git
git -C landslides-trigrs rev-parse HEAD
# Expect 969c409f7bc63e7616ebe34dfd5f6afd967f2d2d
(cd landslides-trigrs/src/TRIGRS && patch -p0 --batch --fuzz=0 -i ../../../release_package/release.patch)
```

`release.patch` differs from the study's registered patch only in explanatory
comments. `logs/patch_verification.log` records zero-fuzz application and the
per-file comment-only diff against the registered release sources: seven files
have comment differences. The correction patch changes ten files relative to
native v2.1.0. The patch was applied to an unmodified copy of the pinned upstream tag. `logs/output_comparison.log` records a
hash-verified registered-release comparison on the long dry-lead regression
cell: every parsed numerical output column agrees exactly.

## Build serial and MPI

Install Python 3, GNU Fortran, and a GNU-Fortran-compatible
[MPICH](https://www.mpich.org/) installation for MPI. `build.py` uses only the
Python standard library, the registered source order and these flags:
`-O2 -std=legacy -fallow-argument-mismatch -ffree-line-length-none`.
Build directories are separate from the sources:

```sh
python3 release_package/build.py landslides-trigrs/src/TRIGRS build_serial
python3 release_package/build.py landslides-trigrs/src/TRIGRS build_mpi --mode mpi --compiler mpif90
```

Both create `trg` in the requested directory. Ensure `mpif90` is the MPICH
wrapper for gfortran (`mpif90 -show`). Use a fresh run directory containing
`tr_in.txt` and its referenced inputs; run the serial executable there, or
use `mpiexec -n 2 build_mpi/trg` for the MPI executable. The upstream
saturated single-output path retains its serial fallback; use serial for
that path. Detailed MPI profile collection is an upstream limitation; this
package's numerical regression is a serial test. Both corrected serial and
MPICH builds were compiled successfully; compiler logs are in `logs/`.

## Run the standalone regression

Install numpy into your Python 3 environment, then:

```sh
python3 release_package/regression/run_regression.py build_serial/trg --workdir regression_run
```

Use a new `--workdir` each time. The runner creates one-cell grids, rainfall
files and `tr_in.txt`, runs the given executable, and checks all requested
times and all eleven depth nodes in the native detailed pressure/FS listing.
It reports the maximum capped-head error (metres) and FS error per case.
Both must be strictly below `1e-4` for PASS. Exit status is 0 for all PASS,
1 for a numerical FAIL, and 2 for a run/input/output error.

The six cases are: an 80-h uniform storm with maximum age about 8; a 6-h
start burst followed by 74 dry hours; an 80-h dry lead and 2-h burst followed
by 3 dry hours with maximum age about 85; two 6-h bursts; an infinite-depth
start burst; and a single-output `savage` case. These small saturated cases
check the finite-depth and half-space serial kernels; they do not establish
independent validation of the coupled unsaturated routines or MPI numerics.

The independent reference uses a mean/storage state and cosine eigenmode
recurrence for the no-flow finite layer, and the Iverson half-space step
response using Python's `math.erfc`. It does not import the project harness
or call a reference executable. It uses the serialized inputs' single
precision, computes the reference in double precision, applies the native
pressure cap and saturated FS convention, and checks mode resolution at
1024 and 2048 modes. Five-significant-figure listing rounding is included
in the error; full-precision instrumentation is unnecessary.

`logs/regression_release.log`: all six cases PASS; largest head and FS
errors are approximately `5e-5`, consistent with listing precision.
`logs/regression_native.log`: the clean, unpatched v2.1.0 serial control
FAILS the five finite-depth cases and PASSES the half-space case.

## Pre-run exposure check

Use `screen/screen.py` before running TRIGRS. It reads the v2.1.0 input file
and referenced grids. Rainfall may come from grids or the nonnegative `cri`
constant overrides. Paths resolve relative to the input file; use
`--base-dir` if the original TRIGRS working directory differs. No solver is
run. Nonnegative `zmax`, `depth` and `rizero` constants replace their property
grids, so those filenames need not exist. A single property zone is assigned
as zone 1 without reading the zone grid. The script reports each property
zone, the slope range actually
present, active/applicable cell counts and the maximum forcing-period age
at the requested output times. Nodata cells are excluded. A given slope
range can replace the observed slopes:

```sh
python3 release_package/screen/screen.py run/tr_in.txt
python3 release_package/screen/screen.py run/tr_in.txt --slopes 10 40
python3 release_package/screen/screen.py landslides-trigrs/tr_in.txt --base-dir landslides-trigrs
```

For a period starting at time `a` and requested output time `t`, age is
`A = D(t-a)/L²`, with `D = D0/cos²(slope)`, vertical soil depth `L` in metres,
diffusivity in m²/s and times in seconds. Count every period, including the
first and periods with zero rain. A slope range is evaluated at its two
endpoints for every cell; ages increase with slope within 0–90 degrees,
so the reported flags conservatively cover the range.

**S21 (late-time exposure).** A finite-depth saturated run is EXPOSED if any
period has age greater than 5 at any requested output time. This also
applies to the unsaturated model's saturated fallback when the stored
water-table depth is at most `1/alpha`. Infinite-depth runs and coupled
unsaturated runs outside that fallback are marked N/A. The screen identifies
entry into the finite-layer late-time branch; it does not assert that the
capped head or FS error will be large.

**S21-D1 (older-period/dry-lead diagnostic).** Report EXPOSED-D1 when an
applicable run has a period aged at least 72 at a requested output time
and a younger nonzero-rain period has already started and is being convolved.
The first old period can be dry: removing it from the check would miss the
dry-lead trigger. Age 72 alone is insufficient. This is a diagnostic for
possible premature history-loop termination, rather than an exact machine
underflow threshold or a requirement that the old period itself be dry.
These are the two v2.1.0 screens in `protocol/scope_protocol_v3.md`.

The tutorial and two harness-generated inputs were checked against the
frozen harness screens. See `logs/screen_tutorial.log` and
`logs/screen_harness.log`. The harness checks cover maximum ages about 1.6
(UNEXPOSED, NO-D1) and 85 (EXPOSED, EXPOSED-D1). Optional-grid tests also cover
absent unused filenames and missing required-grid controls:

```sh
python3 release_package/screen/test_screen.py
```

All eight checks pass; see `logs/screen_optional_grids.log`.

## Licence and archive

This patch, the Python scripts and package documentation are dedicated to
the public domain under [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/),
consistent with the upstream USGS code's public-domain status. See
`LICENSE.md` and the [upstream licence](https://github.com/usgs/landslides-trigrs/blob/v2.1.0/LICENSE.md);
the latter separately identifies `calerf.f` and `derfc.f` as Netlib files.
Retain upstream notices when redistributing the upstream distribution.

The full study archive (inputs, native and corrected outputs, protocols, figure
scripts and this package) is at Zenodo: https://doi.org/10.5281/zenodo.23534840.
