# Actions execution policy — cteasFOAM (2026-10-23)

GitHub still automatically runs `repository-admission` and `scientific-admission` on every PR, with the Ruleset enforcing these checks. The Ubuntu Docker full verified at exact PR SHA is required for source/test PRs. `workflow_dispatch` remains available for the original official/scientific reproducibility workflows; their jobs and tests were not removed.

`docker-dev-image-probe.yml` (30-minute timeout) is dispatch-only, because heavy OpenFOAM compilation is performed locally in Ubuntu Docker at exact PR SHA; the original GitHub job is preserved for explicit independent cross-platform reproductions. This does not imply currently unimplemented HT/HTO, soil PDE or coupling scientific evidence passed.
