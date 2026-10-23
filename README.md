# cteasFOAM

OpenFOAM v2606-based research software for coupled atmospheric and soil tritium transport.

**M0 infrastructure and staged M1A prescribed-flow scalar regressions have
passed their scoped gates. M1B-1f laminar pressure–velocity canonical tests
through three-dimensional warped ABC meshes have passed exact-HEAD v2606 CI.
M1B-1g distorted three-dimensional nonperiodic fixed-pressure duct is an
experimental acceptance candidate. Validated atmospheric URANS, HT/HTO,
soil physics and air–soil coupling are NOT IMPLEMENTED / NOT VALIDATED.**

## Intended architecture

- `src/atmosphere/`: OpenFOAM executable with prescribed-U passive scalar,
  special Couette and experimentally pressure-corrected **laminar** momentum
  modes. URANS closure and HT/HTO activity PDEs are **future milestones**.
- `src/soil/`: independent OpenFOAM 3-D porous-media and HT/HTO solver
  (M0 placeholder only).
- `src/coupling/`: nonmatching 2-D air-soil conservative exchange
  (not implemented).
- `src/models/`, `cases/`, `tests/`, `docker/`, `docs/`: units, case
  contracts, benchmarks and provenance.

## Official OpenFOAM environment

Ubuntu 22.04 + official signed OpenCFD/Keysight OpenFOAM v2606 APT packages,
no third-party OpenFOAM images and no source rebuild. Although v2606 was
officially released 2026-06-26, its official Jammy Debian pool packages are
still named `2606.0~rc2-1`; this exact tested version is intentionally pinned.

On a **native x86-64 Ubuntu Docker host**, execute M0 toolchain checks:

```bash
bash scripts/probe_official_jammy.sh
```

The M1A closed-domain analytic diffusion regression is additionally executed
by GitHub Actions inside the same image. See
`docs/m1a-passive-scalar-acceptance.md` for equation, tolerance and limitations.

## Validation

Compilation is not a physical test. Gate every increment separately with
analytic or conservation benchmarks, field checks, then MPI/scaling and
independent observations. No proprietary data in this public repository.
