# cteasFOAM

OpenFOAM v2606-based research software for coupled atmospheric and soil tritium transport.

**M0 infrastructure and staged M1A prescribed-flow scalar regressions have
passed their scoped gates. M1B-1f laminar pressure–velocity canonical tests
through three-dimensional warped ABC meshes have passed exact-HEAD v2606 CI.
M1B-1g distorted three-dimensional nonperiodic fixed-pressure duct is an
experimental acceptance candidate. M1C-0 adds an **experimental** HT/HTO
activity reaction/transport branch with a closed-box numerical gate. M1C-1
adds a separately gated warped-3D periodic HT/HTO advection–diffusion–reaction
benchmark and finite-volume face-flux activity budget. M1C-2/3/4 add
scope-limited MPI and signed pressure-open activity checks, including real
2-rank x-decomposition and zero-flow controls; none is an independent
scientific or environmental validation. Atmospheric URANS, complete HT/HTO
physics, soil physics and air–soil coupling remain NOT IMPLEMENTED / NOT VALIDATED.**

## Intended architecture

- `src/atmosphere/`: OpenFOAM executable with prescribed-U passive scalar,
  special Couette and experimentally pressure-corrected **laminar** momentum
  modes. Optional M1C HT/HTO advection–diffusion–reaction prototype has
  **restricted numerical** verification in closed boxes, periodic warped
  3-D ducts and open-pressure 3-D ducts, but no independent field validation.
  URANS closure, realistic HT/HTO parameters and deposition are future milestones.
  Optional `auditActivityFlux true` reports the signed HT/HTO matrix
  face flux (Bq/s) and per-species/combined discrete activity residuals;
  M1C-1 verifies a periodic domain. M1C-3 adds a **separately gated**
  nonperiodic pressure-open HT/HTO inflow/outflow and actual finite-volume
  signed Bq budget benchmark, including flow reversal and zero-flow
  negative control. M1C-2 verifies periodic HT/HTO MPI and M1C-4 additionally
  checks serial/MPI pressure-open species fields, signed fluxes and budgets
  for all three boundary modes. M1C-5 introduces optional piecewise-constant
  total HT/HTO release rates (`sourceHT_Bq_s`, `sourceHTO_Bq_s`) integrated
  over the exact overlap of each Euler time step and the dimensioned
  [`sourceStart`, `sourceEnd`) window. The current rate is uniformly distributed
  over the whole actual domain volume; it is **not** a spatially resolved
  stack or accident-release source. See `docs/m1c-5-timed-emission-acceptance.md`
  and `docs/m1c-5-mpi-emission-addendum.md`. These are numerical research
  prototypes, not independent atmospheric/field validation.
- `src/soil/`: independent OpenFOAM 3-D porous-media and HT/HTO solver
  (M0 placeholder only).
- `src/coupling/`: nonmatching 2-D air-soil conservative exchange
  (not implemented).
- `src/models/`, `cases/`, `tests/`, `docker/`, `docs/`: units, case
  contracts, benchmarks and provenance.


## Managed Ubuntu development and merge gates (2026-10-23)

For the shared development host, first read `/home/lhy/dev/AGENTS.md`, `/home/lhy/dev/README.md` and this repository's `AGENTS.md`. Work from the actual GitHub PR HEAD in a clean isolated branch/worktree. Ubuntu Docker runs the primary software tests and exact-SHA numerical checks; GitHub provides the mandatory lightweight `repository-admission` / `scientific-admission` PR checks and selected additional CI. This host operates **Mihomo at `127.0.0.1:7897`**, with no MacBook proxy forwarding. These `/home/lhy/dev` commands are host-specific rather than prerequisites for independent clones.

The cteasFOAM verification environment reuses the pinned local OpenFOAM v2606 image; `local-ci-build` currently supports only CTEAS-next and cufotri, not cteasFOAM.

```bash
/home/lhy/dev/bin/local-ci smoke cteasFOAM --wait
/home/lhy/dev/bin/local-ci full cteasFOAM --wait
/home/lhy/dev/bin/local-ci scientific cteasFOAM --wait
```

The `local-ci` entry points currently require the repository's fixed `/home/lhy/dev/cteasFOAM/worktrees/local-docker-dev` checkout to contain the intended clean HEAD, and image/inputs/results must refer to that SHA. `full` is only the selected **software** suite; `scientific` can legitimately report `PARTIAL`/`NOT_RUN`. A source/test PR needs a verified `ubuntu-docker-full` commit status for its live SHA; the project also retains any applicable original-model, analytical, scientific or independent-review gates.

The `main` Ruleset requires PR review through `repository-admission` and `scientific-admission`. When every applicable check passes and any required independent scientific acceptance is documented, the agent adds `merge-verified` and GitHub's event-driven native auto-merge processes the PR. Draft, failed, pending, missing-data, stale-SHA or scientifically unaccepted changes remain blocked. Heavy scientific workflows can be agent-dispatched without manual GitHub UI interaction; they are **not** automatically PASS when skipped.

## Official OpenFOAM environment

Ubuntu 22.04 + official signed OpenCFD/Keysight OpenFOAM v2606 APT packages,
no third-party OpenFOAM images and no source rebuild. Although v2606 was
officially released 2026-06-26, its official Jammy Debian pool packages are
still named `2606.0~rc2-1`; this exact tested version is intentionally pinned.

On a **native x86-64 Ubuntu Docker host**, execute M0 toolchain checks:

```bash
bash scripts/probe_official_jammy.sh
```

The M1A closed-domain analytic diffusion regression belongs to the Ubuntu Docker smoke/full verification profile; the original GitHub OpenFOAM probe is retained for agent-initiated `workflow_dispatch` reproduction rather than automatic PR builds. See `docs/m1a-passive-scalar-acceptance.md` for equation, tolerance and limitations.

## Validation

Compilation is not a physical test. Gate every increment separately with
analytic or conservation benchmarks, field checks, then MPI/scaling and
independent observations. No proprietary data in this public repository.
