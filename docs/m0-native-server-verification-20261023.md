# M0 native Ubuntu server verification — 2026-10-23

Result: **PASS for the M0 development environment and placeholder build only**.
No atmospheric, soil, interface or parallel numerical model is implemented or validated.

## Tested environment

| Item | Observed value |
| --- | --- |
| Native server | Ubuntu 22.04.5 LTS, x86_64 |
| Container | Official Ubuntu 22.04, linux/amd64 |
| Docker Engine | 29.1.3 |
| OpenFOAM environment | v2606, `linux64GccDPInt32Opt`, `SYSTEMOPENMPI` |
| Official OpenFOAM APT packages | `2606.0~rc2-1` (literal upstream package version) |
| Compiler | GCC/G++ 11.4.0 |
| MPI | Ubuntu Open MPI 4.1.2 (`4.1.2-2ubuntu1`) |
| Ubuntu manifest | `sha256:0a9887931b0e38af2e3487ed39e29436b0535c4f3f52e57745baaf79d1898956` |
| Image content config | `sha256:522281ac6e8f2f70c74f15c515114d5861e9ed325d242292ab4f40756e0b6170` |

OpenFOAM comes exclusively from the signed official `dl.openfoam.com` jammy
APT feed, including its official SourceForge redirects. No third-party
OpenFOAM image, upstream core compilation, or Mac CPU emulation was used.
The Ubuntu manifest and OpenFOAM package versions are pinned in the Dockerfile;
the remaining Ubuntu dependencies are not a complete distribution snapshot.

## Results

| Check | Result / evidence |
| --- | --- |
| Official-package Docker build | PASS; all required development packages installed |
| Required headers | PASS; `fvCFD.H` and `cyclicAMIPolyPatch.H` present |
| Atmosphere placeholder | PASS; fresh wclean/wmake, dynamic linkage and execution |
| Soil placeholder | PASS; fresh wclean/wmake, dynamic linkage and execution |
| FVM / meshTools linkage | PASS; ldd resolves `libfiniteVolume.so`, `libmeshTools.so`, `libOpenFOAM.so`, system OpenMPI Pstream/MPI |
| MPI compiler/runtime | PASS; independent two-rank MPI_Allreduce gives `ranks=2 allreduce=3` |
| Broken MPI version command | PASS; injected exit 42 is preserved and rejected |
| wmake error with exit zero | PASS; injected preprocessing diagnostic is rejected with exit 14 |
| Normal compiler diagnostics | PASS; no compiler errors or missing dependency diagnostics before the explicitly marked negative tests |
| Source provenance | PASS; every tested file SHA-256 matched the local checkout before adding this report |
| Mac host guard | PASS; actual Mac invocation exited 2 with NOT_RUN before Docker use |
| Shell syntax and Git whitespace | PASS; bash -n and git diff --check |
| Solver MPI numerics / scaling | NOT_RUN; not implemented |
| Atmospheric / soil PDE, coupling conservation, field validation | NOT_RUN; not implemented |
| GitHub Actions | Old cc0585cd run 39811328610 failed with NO_PUBKEY 3386DD74948D208F; its key-format bug is fixed here. New-commit CI is separate from native-server evidence |

The negative regressions deliberately print injected errors after the
`EXPECTED FAILURE REGRESSIONS` marker. Their expected rejection statuses are
42 and 14; the complete acceptance probe exits **0**.

## Reproduced issues and fixes

1. The official `.gpg` URL serves an ASCII-armored key. Dearmor it before writing
   a binary `.gpg` keyring; retain `signed-by` scoping and normal APT signatures.
   The subsequently completed old-commit Actions run 39811328610 independently
   reproduced `NO_PUBKEY 3386DD74948D208F` with the old key-handling code.
2. Both Make/options already had the meshTools include/library fix from
   `cc0585cd`. Native compilation confirmed that dependency. Their new `#`
   comments nevertheless caused C-preprocessor errors; change them to C-style
   comments. wmake returned zero despite those diagnostics, so the smoke gate
   now rejects errors as well as nonzero process statuses.
3. MPI version failure was hidden by command substitution inside echo. Invoke
   version commands directly, then compile/run an independent two-rank test.
4. Docker pulls and Buildx client token requests have separate proxy needs.
   The optional build proxy now covers the client and build RUN instructions;
   the daemon pull proxy is configured separately on the server.
5. Keep unique per-attempt logs, source patches and hashes; exclude project
   files from Docker build context and compile only temporary container copies.

## Evidence retained outside Git

Canonical pre-commit run: `official-jammy-20261022T164423Z-bO009m`
(2026-10-23 00:44:23 Asia/Shanghai).

- Full log: `reports/server-before-commit/official-jammy-20261022T164423Z-bO009m/probe.txt`
- Log SHA-256: `d6f1f67903a79f4e2ce7da9e84bc04398ab08bfee59b63a7f6ed8afda0751b2a`
- Source manifest SHA-256: `f5fc2798ae9642dd664d7bb13a60293c7739b078b7ac4dc485cdc7d88c42102f`
- Baseline commit: `cc0585cdafecefa41de65c0b2e2e8bd7fe7be843`; the run retains
  `source.patch` and `source-sha256.txt` to identify all subsequent changes.
- Historical Actions failures: `reports/history/actions-39806603288.txt` and
  `reports/history/actions-39811328610.txt`.
- Earlier server attempts are retained: registry-token timeout (exit 1),
  deliberate cancellation to apply version pins (exit 130), and the pre-fix
  wmake run that returned zero despite preprocessing errors. That pre-fix run
  is **not** accepted as a clean build.

After committing, repeat the same probe on a clean checkout of that commit;
its log identifies the exact commit without a working-tree patch. Retain this
additional evidence in `reports/server-final-commit/` locally. Report the tested
commit SHA in the delivery message rather than making a self-referential commit.

## Reproduction and limits

On a native x86-64 Linux host with Docker Engine and Buildx:

```bash
bash scripts/probe_official_jammy.sh
```

See `docker/README.md` for optional explicit proxy settings and the separate
server/daemon configuration. For this session Docker was installed and started
with the owner's permission. The server's existing Mac reverse-tunnel proxy was
used for downloads; its Docker service proxy override is under `/run` and will
not survive reboot. No permanent docker-group membership was added.

This is M0 infrastructure evidence. The executables still print
`NOT IMPLEMENTED`. The PR remains unmerged; scientific verification and any
formal CI merge gate remain separate decisions.

## Files changed

- `.dockerignore`, `.github/workflows/docker-dev-image-probe.yml`
- `docker/Dockerfile.official-jammy`, `scripts/probe_official_jammy.sh`
- `src/atmosphere/Make/options`, `src/soil/Make/options`
- `tests/smoke/verify_image.sh`, `tests/smoke/mpi_probe.C`, `tests/smoke/verify_failure_paths.sh`
- `README.md`, `docker/README.md`, `docs/docker-image-audit.md`, this report
