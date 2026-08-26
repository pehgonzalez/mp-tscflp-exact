# MP-TSCFLP, exact algorithms and benchmark results

Code, instances and raw results for the computational study of the multiproduct
two-stage capacitated facility location problem (MP-TSCFLP). The repository
contains everything needed to rebuild the solver, rerun the campaigns and audit
the reported numbers.

## Contents

| Path | Content |
| --- | --- |
| `src/` | C++ sources of the solver. One binary, `mptscfl`, drives the three engines. |
| `tests/` | Solver-free core tests. |
| `data/instances/` | The 100 benchmark instances. `FORMAT_original.txt` documents the file format. |
| `data/bks_reference.csv` | Published best-known values used as reference. |
| `scripts/` | Campaign drivers for Windows PowerShell. `run_campanha_revisao.ps1` runs the full revision campaign, `run_extensao_rodada3.ps1` the certification and factorial extension with a progress dashboard. |
| `python/` | Table, figure and statistics generators, reference implementations in gurobipy, and the verification tools. |
| `results/results.csv` | The complete campaign record, 1,620 runs, one row per run. |
| `results/certdump/` | The recorded final masters of the six audited closures. |
| `results/rational_report_v2.txt` | Output of the exact-arithmetic audit over those six masters. |

## Build

Requires CMake 3.16+, a C++20 compiler and Gurobi. The study ran on Windows
with Gurobi 13.0.2. Set `GUROBI_HOME` or pass `-DGUROBI_HOME=...`.

Windows, the platform of the study:

```powershell
cmake -B build -DGUROBI_HOME=C:/gurobi1300/win64
cmake --build build --config Release
ctest --test-dir build -C Release
```

Linux:

```sh
cmake -B build -DGUROBI_HOME=/opt/gurobi1300/linux64
cmake --build build
ctest --test-dir build
```

Without Gurobi only the core tests build. On Windows the binary lands at
`build\Release\mptscfl.exe`, on Linux at `build/mptscfl`.

## Run

```
mptscfl <instance> <method> <time_limit_s> [mode] [seed] [threads] [corepoint]
```

| Argument | Meaning |
| --- | --- |
| `instance` | Path to an instance file from `data/instances/`. |
| `method` | The engine. `0` solves the compact model directly with the MIP solver, `1` is branch-and-Benders-cut (BBC), `2` is the Lagrangian-guided variant (L-BBC). |
| `time_limit_s` | Time budget in seconds for the run. |
| `mode` | Configuration label, default `two-steps`. The table below lists every mode. |
| `seed` | Gurobi `Seed`, default 0. Every recorded run of the study uses seed 0. |
| `threads` | Gurobi `Threads`, default 0 meaning automatic. The study used 16. |
| `corepoint` | `1` turns on core-point cuts (Magnanti and Wong in the Papadakos variant) for the Benders engines, default 0. Used together with the label `exact-pap`. |

Every run writes a full solver log to
`logs/<instance>_m<method>_<mode>_s<seed>_<timestamp>.gurobi.log` and appends
one row to `logs/results.csv`. The row includes `verified_cost` and
`verified_ok`, a solver-free re-routing of the returned opening decisions that
recomputes the cost independently before it is reported.

### Modes

The mode is a label. `two-steps` is the default. Any label starting with
`exact` runs the exact configuration, and the suffix selects an ablation or a
protocol variant. The `runs` column counts the rows each label has in
`results/results.csv`.

| Mode | Engine | What it does | Runs |
| --- | --- | --- | --- |
| `two-steps` | 0, 1, 2 | Runs the two-steps decomposition heuristic first and hands its solution to the chosen engine as a starting point and cutoff. Default mode, not part of the recorded campaign. | 0 |
| `exact` | 0, 1, 2 | The canonical configuration of the study. | 790 |
| `exact-var1` `-var2` `-var3` | 0, 1, 2 | Identical to `exact`. The labels keep repetition runs apart in the record. | 51 |
| `exact-pap` | 1, 2 | The core-point ablation. Run with the last argument set to 1. | 133 |
| `exact-agg` | 1 | One aggregated optimality cut instead of one cut per product. | 25 |
| `exact-noldb` | 2 | Withholds the Lagrangian bound inequality from the master. The repaired start and the cutoff stay. | 77 |
| `exact-nowarm` | 2 | The subgradient phase starts at zero multipliers instead of the linear programming duals. | 50 |
| `exact-nostart` | 2 | Keeps the bound inequality and the cutoff, withholds the repaired starting solution. | 50 |
| `exact-oneeval` | 2 | Reduces the Lagrangian phase to a single evaluation of the dual function at the linear programming duals, with no multiplier update. | 100 |
| `exact-fl` | 0 | The compact model plus the two aggregate feasibility rows per product. | 100 |
| `exact-mipguided` | 0 | The compact model with the recorded Lagrangian bound as an objective floor, the repaired solution as a MIP start, and the cutoff. | 100 |
| `exact-flguided` | 0 | The two previous rows combined, feasibility rows together with the guidance. | 100 |
| `exact-mauri` | 0 | Solver default tolerances and no proof mode, the protocol of the study that introduced the benchmark. | 12 |
| `exact-cert` `-cert2` `-cert3` | 1, 2 | Labels of the certificate-dump reruns. Run with the environment variable `MPTSCFL_CERT_DUMP` set, see below. | 32 |
| `exact-repaircb` | 0 | The compact model with the incumbent-repair callback, which re-routes every candidate and injects the cheaper repaired solution. Available in the code, no recorded runs. | 0 |
| `exact-flgcb` | 0 | `exact-flguided` plus that same callback. Available in the code, no recorded runs. | 0 |
| `exact-thetaint` | 1, 2 | Declares the value variables integer in the master. Requires the source patch described in `src/main.cpp` and refuses to run without it. No recorded runs. | 0 |

### Environment variables

| Variable | Effect |
| --- | --- |
| `MPTSCFL_GRB_PARAMS` | Extra Gurobi parameters injected into every model the run builds, as `Name=Value` pairs. |
| `MPTSCFL_CERT_DUMP` | A directory path. The Benders engines record the final master there as a `.certdump` file, cuts, incumbent and bounds, which is the input of the exact-arithmetic audit. |

### Examples

All examples use PowerShell paths. On Linux, replace the binary path with
`build/mptscfl` and the backslashes with slashes.

The canonical run of the study, L-BBC at the thirty-minute budget:

```powershell
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 2 1800 exact 0 16 0
```

The same instance under the other two engines, plain BBC and the direct solve
of the compact model:

```powershell
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 1 1800 exact 0 16 0
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 0 1800 exact 0 16 0
```

The default mode, a two-steps heuristic solution handed to the engine as start
and cutoff:

```powershell
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 1 1800 two-steps 0 16 0
```

The factorial arms of the study, each removing one component of the guidance:

```powershell
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 2 1800 exact-noldb 0 16 0
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 2 1800 exact-nowarm 0 16 0
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 2 1800 exact-nostart 0 16 0
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 2 1800 exact-oneeval 0 16 0
```

The core-point ablation, label plus flag:

```powershell
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 1 1800 exact-pap 0 16 1
```

The single-cut master ablation:

```powershell
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 1 1800 exact-agg 0 16 0
```

The compact-model variants, feasibility rows, guidance, and both at once:

```powershell
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 0 1800 exact-fl 0 16 0
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 0 1800 exact-mipguided 0 16 0
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 0 1800 exact-flguided 0 16 0
```

The rerun under the protocol of the study that introduced the benchmark, one
hour and solver default tolerances:

```powershell
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 0 3600 exact-mauri 0 16 0
```

A certificate-dump rerun, which records the final master for the audit:

```powershell
$env:MPTSCFL_CERT_DUMP = "results\certdump"
build\Release\mptscfl.exe data\instances\PSC1-C1-50-5.txt 1 14400 exact-cert3 0 16 0
```

The campaign drivers resume from `results.csv` and skip runs already recorded:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_campanha_revisao.ps1 -Mode Report
powershell -ExecutionPolicy Bypass -File scripts\run_campanha_revisao.ps1 -Mode Apply
```

## Verify the reported numbers

Every table and figure of the paper regenerates from `results/results.csv`.

```sh
python3 python/make_tables.py results/results.csv data/bks_reference.csv
```

The exact-arithmetic audit of the six recorded masters runs with the standard
library only and needs no solver.

```sh
python3 python/rational_check.py results/certdump/PSC1-C1-50-5_20260819-222709.certdump --instance data/instances/PSC1-C1-50-5.txt
python3 python/rational_check.py --selftest
```

## License

MIT, see `LICENSE`. The benchmark instances remain the intellectual property of
their original authors.
