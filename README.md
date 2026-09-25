# Dynein Transport Simulation

This repository contains the C simulation code used in the associated article, **“Modelling nano-particle-PEG-NLS complexes for nucleus-targeted drug delivery: revisiting the multi-dynein nano-cargo transport on microtubules”**, by Tal Yagev, Itay Fayer, Itay Adar, Gal Halbi, Anne Bernheim-Groswasser, and Rony Granek.

## Scientific Overview

The simulation models the transport of a spherical nanoparticle along a microtubule by multiple dynein motors connected to the particle through flexible polymer linkers. It was developed to study PEG–NLS-decorated nanoparticles that recruit dynein for transport toward the nuclear envelope.

The model extends our previous multi-dynein transport framework by introducing two key features: **semi-flexible dynein motors** and **steric exclusion along the complete motor-stepping path**. Motor binding, unbinding, and stepping are simulated stochastically, while the nanoparticle, linkers, and motors mechanically relax after each event through energy minimization.

The simulations are used to study how **dynein loading, PEG linker length, and nanoparticle properties affect transport velocity, run distance, and run time**. The model reproduces the experimentally observed dependence of nanoparticle motility on motor number and can therefore be used to explore design rules for dynein-powered, nucleus-directed nanoparticle delivery. Under the conditions studied in the paper, the simulations predict an optimal PEG contour length of approximately **40 nm** for maximizing run distance and processivity.

The simulation writes its results to `dynein_simulation_output.txt`. Full details of the model, parameters, assumptions, and analysis are provided in the manuscript and Supporting Information.

## Files

| Path | Description |
| --- | --- |
| `src/dynein_transport_simulation.c` | Main simulation source code. |
| `dynein_simulation_input.txt` | Easy-to-edit input file for the simulation parameters. |
| `scripts/` | Convenience scripts for building and running the simulation. |
| `docs/` | Setup notes. |

## Input Parameters

Edit `dynein_simulation_input.txt` before running the simulation. **The file contains one value per parameter, each parameter in a separate line, in the following order:**

1. Number of grafted dynein motors on the nanoparticle
2. Number of nanoparticle simulation runs
3. Nanoparticle radius in meters
4. Dynein persistence length in meters
5. Polymer length in meters

Lines starting with `#` are comments and are ignored by the simulation.

All three length parameters must be entered in meters, not nanometers or micrometers. For example, 20 nm is `2.0e-8` meters, and 2 micrometers is `2.0e-6` meters. The number of grafted motors is distinct from the number simultaneously bound to the microtubule, which changes during a simulated run.

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

This creates:

```text
bin/dynein_transport_simulation
```

### Windows PowerShell

If `gcc` is available in your PATH, run:

```powershell
.\scripts\build_all.ps1
```

This creates:

```text
bin\dynein_transport_simulation.exe
```

For Visual Studio, create a C console project and add `src/dynein_transport_simulation.c`. See `docs/windows_pthreads.md` for pthreads setup notes.

## Run

Make sure `dynein_simulation_input.txt` is in the repository root, then run from that directory:

```sh
./bin/dynein_transport_simulation
```

On Windows PowerShell:

```powershell
.\bin\dynein_transport_simulation.exe
```

## Output and Interpretation

The simulation writes:

```text
dynein_simulation_output.txt
```

The scientific analysis in the associated manuscript concerns nanoparticle transport statistics, including longitudinal velocity, run distance, run time, the number of microtubule-bound motors, and rotational motion. These observables are defined in the manuscript; their calculation may involve analysis beyond the raw simulation output.

Keep each output file together with the corresponding input file and code version so that the simulation conditions remain identifiable.

## Reproducibility Notes

- The simulation uses stochastic sampling, so repeated runs may differ unless random seeding is made fully explicit.
- Runtime depends strongly on the number of simulation runs and motors.
- The five input-file parameters do not constitute the entire model specification. Consult the source code, manuscript, and supporting information for other model constants and assumptions.
- Compiled binaries, object files, and temporary output files are intentionally excluded from version control.

## Citation

If you use this software, please cite:

> Tal Yagev, Itay Fayer, Itay Adar, Gal Halbi, Anne Bernheim-Groswasser, and Rony Granek (2026). “Modelling nano-particle-PEG-NLS complexes for nucleus-targeted drug delivery: revisiting the multi-dynein nano-cargo transport on microtubules.” *Drug Delivery*. [https://doi.org/10.1080/10717544.2026.2738276](https://doi.org/10.1080/10717544.2026.2738276)

Citation metadata are also available in [`CITATION.cff`](CITATION.cff).

Repository: [Dynein Transport Simulation](https://github.com/tyagev-dotcom/dynein-transport-simulation).

## License

The simulation code is released under the BSD 3-Clause License. See `LICENSE.md`.

The software license does not apply to the associated manuscript or supporting information, which retain their applicable copyright and sharing terms.
