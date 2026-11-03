# Minute-scale GPU benchmark v1

This fixed workload measures a complete MBS Mueller calculation after at least
120 seconds of excluded warmup. It uses a hexagonal prism at ka=10, n=1.3116,
2048 deterministic Sobol/Haar orientations, **2048 azimuth samples**, and 181
scattering angles with all 16 matrix elements. Angular work is 32 times larger
than the short `hardware_v1` case. It does not change the physics solver.

```bash
python3 scripts/benchmark_gpu_thermal.py run --output runs/thermal-5080
```

The default profile is explicit `gpu-fp32-reduced`; `--profile gpu-fp32` and
`--profile gpu-fp64` use the precise phase-storage profiles. The separately
available double-single profile has not been scored in this suite. Install
the desired binaries with [the RTX 50 installer](../../docs/INSTALL_RTX50.md).
Use `--prefix` to select an installation. The runner records its actual native
binary, compiler, source and toolkit hashes, including when that installation
was built from an earlier revision than the runner.

The runner fixes eight CPU threads, performs at least three fresh measured
processes spanning at least 180 seconds, and validates every full curve against
[the frozen FP64 reference](reference/fp64.json). Regenerate a candidate
reference explicitly with `reference --output NEW_DIRECTORY`; it does not
replace the published reference automatically.

A measured run is eligible only with native CUDA kernel time/wall time >=90%,
median GPU load >=90%, and passing matrix/work-count gates. CPU tracing, process
startup and output remain in wall time. The current runner also samples compute
PIDs throughout each invocation and rejects foreign or unidentifiable GPU jobs.
Run on an idle GPU; do not compare shared-GPU measurements. It does not stop
other jobs. Raw per-second load, temperature, power and clock telemetry are
preserved alongside CPU time, RSS, phase timings and exact commands.

```bash
python3 scripts/benchmark_hardware.py compare baseline/results.json candidate/results.json
python3 benchmarks/gpu_thermal_v1/plot_results.py runs/thermal-5080/results.json
```

Compare only equal workload/profile/reference/protocol keys. Source/toolchain
differences are reported rather than interpreted as pure hardware speedups.
These are discrete numerical regression checks, not a certificate of physical
convergence with respect to orientation, angle sampling or PO approximation.

[RTX 5080 measurement and provenance](results/rtx5080_20261103/README.md).
