# Dynein Transport Simulation

This repository contains a C simulation for dynein-driven transport of a rigid nanoparticle along a microtubule.

The code was prepared for manuscript revision and public archival. The repository should remain private until the manuscript team approves public release.

## Files

| Path | Description |
| --- | --- |
| `src/dynein_transport_simulation.c` | Main simulation source code. |
| `dynein_simulation_input.txt` | Easy-to-edit input file for the simulation parameters. |
| `scripts/` | Convenience scripts for building and running the simulation. |
| `docs/` | Setup notes and publication checklist. |

## Input Parameters

Edit `dynein_simulation_input.txt` before running the simulation. The file contains one value per parameter:

1. Number of grafted dynein motors on the nanoparticle
2. Number of nanoparticle simulation runs
3. Nanoparticle radius in meters
4. Dynein persistence length in meters
5. Polymer length in meters

Lines starting with `#` are comments and are ignored by the simulation.

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

Make sure `dynein_simulation_input.txt` is in the repository root, then run:

```sh
./bin/dynein_transport_simulation
```

On Windows PowerShell:

```powershell
.\bin\dynein_transport_simulation.exe
```

The simulation writes:

```text
dynein_simulation_output.txt
```

## Reproducibility Notes

- The simulation uses stochastic sampling, so repeated runs may differ unless random seeding is made fully explicit.
- Runtime depends strongly on the number of simulation runs and motors.
- Compiled binaries, object files, local PDFs, and temporary output files are intentionally excluded from version control.

## Citation

Please cite the associated manuscript when using this code. A `CITATION.cff.template` file is included and should be finalized after the manuscript metadata are confirmed.

## License

No open-source license has been selected yet. Before public release, the manuscript team should choose a license and replace `LICENSE.md` with the final approved license text.
