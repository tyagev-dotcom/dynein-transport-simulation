# Dynein Nanoparticle Transport Simulation

This repository contains C implementations of stochastic simulations for dynein-driven nanoparticle transport along a microtubule.

The code was prepared for manuscript revision and public archival. The repository is currently intended to remain private until the manuscript team approves public release.

## Repository Contents

| Path | Description |
| --- | --- |
| `src/rigid_nanoparticle.c` | Baseline rigid nanoparticle simulation. |
| `src/rigid_nanoparticle_low_persistence.c` | Rigid nanoparticle simulation with modified dynein persistence length. |
| `src/rigid_nanoparticle_flexible_input.c` | Rigid nanoparticle simulation that reads key physical parameters from `Input-NP-flex.txt`. |
| `src/fluid_nanoparticle.c` | Fluid/lipid nanoparticle simulation. |
| `src/fluid_nanoparticle_optimized.c` | Optimized fluid/lipid nanoparticle simulation. |
| `examples/Input-NP-flex.txt` | Example input file for `rigid_nanoparticle_flexible_input.c`. |
| `examples/reference_outputs/` | Small reference output files from prior runs. |
| `reference/original_final_project_code.c` | Original final project script kept for provenance. |
| `scripts/` | Convenience scripts for building and running simulations. |
| `docs/` | Setup notes and publication checklist. |

## Requirements

- C compiler with C11 support
- POSIX threads support
- Math library support

On Linux, this is usually available through `gcc` and `pthread`.

On Windows, use one of:

- MSYS2/MinGW-w64 with `gcc` and `winpthreads`
- Visual Studio plus a pthreads-win32 compatible library

## Build

### Linux or MSYS2

From the repository root:

```sh
make
```

This creates executables in `bin/`.

To build one target:

```sh
make rigid
make rigid-flex
make fluid-opt
```

### Windows PowerShell

If `gcc` is available in your PATH, run:

```powershell
.\scripts\build_all.ps1
```

For Visual Studio, create a C console project and add the desired file from `src/`. See `docs/windows_pthreads.md` for pthreads setup notes.

## Run

Example:

```sh
./bin/rigid_nanoparticle
```

For the flexible-input simulation, copy the example input file to the repository root before running:

```sh
cp examples/Input-NP-flex.txt Input-NP-flex.txt
./bin/rigid_nanoparticle_flexible_input
```

The simulations write output text files in the working directory. The exact output file name is currently defined inside each C source file.

## Simulation Variants

| Executable | Source file | Output file |
| --- | --- | --- |
| `rigid_nanoparticle` | `src/rigid_nanoparticle.c` | `output-NP.txt` |
| `rigid_nanoparticle_low_persistence` | `src/rigid_nanoparticle_low_persistence.c` | `output-NP-LP1.txt` |
| `rigid_nanoparticle_flexible_input` | `src/rigid_nanoparticle_flexible_input.c` | `output-NP-flex.txt` |
| `fluid_nanoparticle` | `src/fluid_nanoparticle.c` | `output-NLP.txt` |
| `fluid_nanoparticle_optimized` | `src/fluid_nanoparticle_optimized.c` | `output-NLP.txt` |

Note: `fluid_nanoparticle.c` and `fluid_nanoparticle_optimized.c` both write to `output-NLP.txt`; run them in separate folders or rename the output in the source before comparing runs.

## Reproducibility Notes

- The simulation uses stochastic sampling, so repeated runs may differ unless random seeding is made fully explicit.
- Some simulations can run for minutes to hours depending on the number of trials and model variant.
- Compiled binaries, object files, local PDFs, and temporary output files are intentionally excluded from version control.

## Citation

Please cite the associated manuscript when using this code. A `CITATION.cff.template` file is included and should be finalized after the manuscript metadata are confirmed.

## License

No open-source license has been selected yet. Before public release, the manuscript team should choose a license and replace `LICENSE.md` with the final approved license text.
