# Agent development instructions

## Owner's current CI and merge decision (2026-10-23)
- Retain automatic GitHub PR admission gates; the heavy OpenFOAM v2606 image/M1A probe is now retained as agent-dispatched `workflow_dispatch`, not an automatic PR build. No owner-initiated manual click is required; prior blanket manual-only Draft PR #21 was closed.
- Continue Ubuntu Docker compile/smoke/full and targeted analytic/conservation checks at the exact candidate SHA; relevant lightweight GitHub PR admission checks then trigger automatically, while heavy original OpenFOAM probes can be explicitly dispatched.
- After every applicable software, analytic, independent validation and GitHub check passes, automatically mark PR Ready and squash merge or enable native auto-merge without a separate owner request. Do not infer HT/HTO, soil physics, coupling or URANS scientific PASS from existing M0/M1 tests. Draft PR #20 remains unmerged until its actual scope is accepted.


Project: RRain001/cteasFOAM. All code is public. Never upload proprietary experimental data, restricted PDFs, personal credentials, or copyrighted upstream OpenFOAM source snapshots without rights review.


## Ubuntu Docker / GitHub admission (current deployment, 2026-10-23)
- On the managed Ubuntu host, read `/home/lhy/dev/AGENTS.md` and `/home/lhy/dev/README.md` **before** this repository guide; in other environments those absolute paths may not exist, so use the equivalent documented container and gate contracts. A new ChatGPT session must explicitly read the files, not assume they are automatically injected.
- Refresh remote Git refs and fetch the exact PR HEAD SHA; inspect clean worktrees and all existing tests. Use an isolated short-lived task branch/worktree, never commit directly to `main`, and never force-clean an accepted worktree. The Ubuntu working copy is under `/home/lhy/dev/cteasFOAM/worktrees/local-docker-dev`; the `local-ci` tool currently binds that path, so check the selected SHA before each run.
- Run real builds, focused smoke, scoped `full`, and admissible scientific checks on Ubuntu in Docker. Host entry points: `/home/lhy/dev/bin/local-ci-build`, `/home/lhy/dev/bin/local-ci`, `/home/lhy/dev/bin/local-ci-report`; follow `/home/lhy/dev/infra/local-ci/README.md` for flags and resource limits. A scoped software pass is not all-option or independent scientific acceptance.
- Ubuntu runs its **own Mihomo** HTTP(S) proxy at `127.0.0.1:7897`. Do not require a MacBook/Clash reverse tunnel. Docker daemon, BuildKit and the host-side buildx CLI need compatible proxy settings. Never commit proxy subscriptions, credentials or restricted data.
- The GitHub `main` Ruleset requires the always-running `repository-admission` and `scientific-admission` PR checks. For applicable source/test PRs, publish `ubuntu-docker-full` success **only** after verified `full` software evidence at the exact live PR SHA via `python3 /home/lhy/dev/bin/local-ci-publish-status cteasFOAM <pr> <full-run-id>`; retrigger admission when necessary. This status does not assert science PASS.
- Agents can dispatch longer workflows at a reviewed candidate SHA without asking the owner to click GitHub Actions. Record the run's actual SHA, inputs, checks and scientific conclusion. For science changes, obtain the project-specific external evidence and a genuinely independent, same-SHA human/non-bot APPROVED review as required; never manufacture a `science-validated` claim.
- Only after **all applicable** checks and any required source-specific scientific evidence and independent review succeed, and the current PR is non-Draft and conflict-free, may the agent add `merge-verified`. The trusted `pull_request_target` workflow enables GitHub native squash auto-merge; changing the SHA invalidates the old authorization. NOT_RUN, FAIL, pending, Draft and scientific-review blockers must stay unmerged. No hourly polling or manual-owner click is required.

## Fixed baseline and boundaries
- OpenFOAM **OpenCFD/Keysight v2606**, Ubuntu 22.04 x86-64 Linux server, Docker-preferred.
- Two **separate executables**, both linked against OpenFOAM v2606 C++ FVM libraries:
  - 3-D atmospheric URANS + HT and HTO transport.
  - 3-D porous soil PDE: HT oxidation, water/HTO transport, gas–liquid exchange.
- Actual 2-D air–soil physical interface with nonmatching meshes, conservative HT/HTO activity exchange and multirate synchronization.
- Plant, food chain and dose are separately validated kinetic models.
- Do not rewrite OpenFOAM FVM/MPI infrastructure or fork upstream until a specific upstream limitation is demonstrated.
- Source of truth for thesis scope and timing remains private `career-plan`; do not copy private documents here.

## Verification discipline
- Mark skeletons as **NOT IMPLEMENTED**. Never imply a solver exists because `wmake` succeeds.
- Progress must be evidenced: compilation -> analytic/conservation tests -> independent field validation -> HPC scaling.
- HT and HTO have independent activity balances; net conversion preserves total tritium activity before radioactive decay and boundary fluxes.
- Positive interface sign is **air to soil**. Transfer signed integrated activity in Bq, not just concentration, on each coupling window.
- On public GitHub, do not import restricted MUST/JAEA/other nonredistributable experimental documents.

## M0: only official OpenFOAM packages
- Build `docker/Dockerfile.official-jammy` from `ubuntu:22.04` and install the precompiled OpenCFD v2606 development packages from the official signed APT repository.
- Verify `wmake`, `g++`, `mpirun`, OpenFOAM headers, and actual compilation of two clearly labeled placeholder applications.
- Do not rebuild the complete OpenFOAM source tree unless a demonstrated issue requires it.
- Do not add external OpenFOAM container distributions as an alternate deployment path without a new explicit decision.
