# Docker strategy

No from-scratch OpenFOAM build for M0.

First inspect official `opencfd/openfoam-dev:2606`. If the tag is missing or the compilation smoke fails, inspect `microfluidica/openfoam:2606`, which documents that its non-slim image can build OpenFOAM applications. **Never use `opencfd/openfoam-run:2606` to compile**: it is runtime-only.

## On an Ubuntu Linux host

```bash
bash scripts/probe_dev_images.sh
```

The script downloads candidate images, checks the toolchain and uses `wmake` to compile and run M0 placeholders in a temporary container directory. Successful smoke tests alone do not prove flow or reactive-transport correctness.

After selecting an image, pin it by digest, keep OpenFOAM versions identical for both solvers, and mount cases/results outside the container writable layer. Run MPI scalability tests separately.
