# MP-TSCFLP, exact algorithms and benchmark results

Code, instances and raw results for the computational study of the multiproduct
two-stage capacitated facility location problem (MP-TSCFLP). The repository
contains everything needed to rebuild the solver, rerun the campaigns and audit
the reported numbers.

## Contents

| Path | Content |
| --- | --- |
| `src/` | C++ sources. Method 0 is the direct MIP solve, 1 branch-and-Benders-cut, 2 the Lagrangian-guided variant. |
| `tests/` | Solver-free core tests. |
| `data/instances/` | The 100 benchmark instances, distributed with permission of their original authors. `FORMAT_original.txt` documents the file format. |
| `data/bks_reference.csv` | Published best-known values used as reference. |
| `scripts/` | Campaign drivers for Windows PowerShell. `run_campanha_revisao.ps1` runs the full revision campaign, `run_extensao_rodada3.ps1` the certification and factorial extension with a progress dashboard. |
| `python/` | Table, figure and statistics generators, reference implementations in gurobipy, and the verification tools. |
| `results/results.csv` | The complete campaign record, 1,620 runs, one row per run. |
| `results/certdump/` | The recorded final masters of the six audited closures. |
| `results/rational_report_v2.txt` | Output of the exact-arithmetic audit over those six masters. |

## Build

Requires CMake 3.16+, a C++20 compiler and Gurobi. Set `GUROBI_HOME` or pass
`-DGUROBI_HOME=...`.

```sh
cmake -B build -DGUROBI_HOME=/opt/gurobi1300/linux64
cmake --build build --config Release
ctest --test-dir build
```

Without Gurobi only the core tests build.

## Run

A single run takes the instance, the method, the time limit in seconds and, as
positional extras, the mode, the seed, the thread count and the core-point flag.
The canonical configuration of the study is mode `exact`, seed 0, core-point off.

```sh
build/mptscfl data/instances/PSC1-C1-50-5.txt 2 1800 exact 0 16 0
```

The campaign drivers resume from `results.csv` and skip runs already recorded.

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
their original authors and are redistributed here with permission.
