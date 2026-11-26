# TITAN RTX resident periodic eigensolver · archived measurement

On 16 October 2026 the four-cell suite ran on a TITAN RTX and dual Xeon Silver 4210 host, using FEM revision `5bd9fd1`. Each iris, substrate-integrated waveguide (SIW), perforated guide and dielectric-rod mesh was assembled once. CPU ARPACK, hybrid GPU and resident GPU solved the same reduced quadratic pencil. Assembly and subsequent field reconstruction are outside the timer.

The suite requested four eigenpairs with ncv = 80 and tolerance 10⁻⁸. Each configuration used one warmup on each backend and three timed repeats with rotated order, yielding 171 timed solves and 57 warmups. All 19 configurations validated. Maximum QEP residual was 2.84 × 10⁻¹³. The Linux/TITAN hardware run is preserved here for historical comparison.

| Case | Target (mm) | Field DOFs | CPU (s) | Hybrid (s) | Resident (s) | CPU / resident |
| --- | --- | --- | --- | --- | --- | --- |
| iris-coarse | 4.0 | 1,591 | 0.2165 | 0.1704 | 0.2079 | 1.04× |
| iris-medium | 3.0 | 2,422 | 0.3699 | 0.1721 | 0.2144 | 1.73× |
| iris-fine | 2.0 | 5,477 | 1.154 | 0.2259 | 0.2604 | 4.43× |
| iris-large | 1.4 | 11,782 | 4.071 | 0.3862 | 0.3656 | 11.1× |
| siw-coarse | 2.0 | 2,145 | 0.2651 | 0.1613 | 0.1984 | 1.34× |
| siw-medium | 1.4 | 3,646 | 0.4955 | 0.1807 | 0.2317 | 2.14× |
| siw-fine | 1.0 | 6,623 | 1.148 | 0.2291 | 0.2626 | 4.37× |
| siw-large | 0.75 | 11,711 | 2.876 | 0.3544 | 0.3745 | 7.68× |
| perforated_waveguide-coarse | 5.0 | 8,373 | 1.982 | 0.3132 | 0.3588 | 5.52× |
| perforated_waveguide-medium | 3.5 | 9,685 | 2.449 | 0.3179 | 0.322 | 7.61× |
| perforated_waveguide-fine | 2.5 | 19,680 | 6.86 | 0.578 | 0.6039 | 11.4× |
| perforated_waveguide-large | 2.0 | 32,161 | 18.45 | 0.7215 | 0.7144 | 25.8× |
| dielectric_rod-coarse | 4.0 | 2,484 | 0.6657 | 0.1785 | 0.2149 | 3.1× |
| dielectric_rod-medium | 3.0 | 3,720 | 1.231 | 0.2088 | 0.2264 | 5.44× |
| dielectric_rod-fine | 2.0 | 7,807 | 5.084 | 0.3331 | 0.3619 | 14.0× |
| dielectric_rod-large | 1.7 | 11,609 | 15.24 | 0.4372 | 0.4116 | 37.0× |
| perforated_waveguide-xl | 1.7 | 45,926 | 34.36 | 1.017 | 1.043 | 32.9× |
| perforated_waveguide-xxl | 1.3 | 83,298 | 119.0 | 1.778 | 1.845 | 64.5× |
| perforated_waveguide-100k | 1.15 | 108,261 | 198.4 | 2.299 | 2.454 | 80.9× |

![Measured TITAN RTX speedups](../../benchmarks/resident_comparison/results/2026-10-16/speedups.png)

[Source report](https://github.com/SolverNotConverging/FEM/blob/main/benchmarks/resident_comparison/results/2026-10-16/REPORT.rst) · [CSV](../../benchmarks/resident_comparison/results/2026-10-16/comparison.csv) · [Audit](../../benchmarks/resident_comparison/results/2026-10-16/audit.json) · [Reproduction instructions](https://github.com/SolverNotConverging/FEM/blob/main/benchmarks/resident_comparison/README.rst)

[Benchmark index](README.md)
