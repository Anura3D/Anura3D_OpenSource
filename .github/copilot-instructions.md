<!--
Guidance for AI coding assistants working on the Anura3D_OpenSource repository.
Keep this short, concrete and focussed on repository-specific patterns and developer workflows.
-->

# Quick instructions for AI coding agents

Be concise and prefer small, local edits with targeted tests. This project is a Fortran-based scientific code (Material Point Method). Key facts you will need immediately:

- Language & build: primarily Fortran (fixed/free form) with Intel Fortran project files under `src/VS/`. The canonical entrypoint is `src/Anura3D.for` (program Anura3D -> calls InitialiseDimension() and Kernel()). Visual Studio solution `src/VS/Anura3D.sln` and `.vfproj` files target Intel Fortran and link to Intel MKL (mkl_dss).

- Major subsystems (high level):
  - Kernel & initialization: `src/Kernel.for`, `src/InitialiseKernel.FOR`, `src/Anura3D.for`.
  - Input/Output & file parsing: `src/FileIO.for`, `src/ReadCalculationData.for`, `src/ReadGeometryData.FOR`, `src/ReadMaterialData.for` (GOM reader).
  - Material models: `src/Soilmodels/` (see `Readme.md` inside that dir). Example models: `UMAT_LinearElasticity.f`, `A3DMohrCoulombStandard.f`.
  - Solvers & third-party: `src/Solver/` (uses Intel MKL; `mkl_dss.f90` is included). `src/VS/*` contains Visual Studio project settings tuned for Intel compilers.

- Platform & dependency notes:
  - The codebase expects Intel Fortran compiler features (conditional compilation flags such as __INTEL_COMPILER, and modules like IFCore/IFPort). Build instructions commonly use Visual Studio solution files on Windows.
  - MKL Pardiso/DSS is used in `src/Solver/Solver.FOR`. Treat MKL calls as external — avoid modifying MKL-specific headers unless necessary.

- Material model extension pattern (important example):
  - External soil models (UMAT-style) live in `src/Soilmodels/`. Follow `src/Soilmodels/Readme.md` when adding models: add module wrapper, implement `ESM_*` subroutine interface, and register the pointer in `src/ReadMaterialData.for` (uncomment assignment like `MatParams(I)%ESM_POINTER => ESM_CamClay`).
  - Example files: `src/Soilmodels/UMAT_LinearElasticity.f` and `src/Soilmodels/A3DLinearElasticity.f` show the expected `module` / `ESM` signature and usage of arrays `ESM_Solid` and `ESM_Statvar_in`.

- Coding conventions & patterns to follow (discoverable):
  - Modules are used per subsystem (e.g., `ModReadMaterialData`, `ModGlobalConstants`). Prefer `use` statements as other files do.
  - Fixed/free-form Fortran mixed across files; preserve original formatting and line length conventions when editing.
  - Error reporting uses `GiveError`, `GiveMessage`, and `Assert` helpers (look in `ErrorHandler.for`/`Feedback.for`). Use these utilities for consistent diagnostics.

- Developer workflows (what an agent might need to suggest):
  - Typical local dev build: open `src/VS/Anura3D.sln` in Visual Studio (Intel Fortran). The `.vfproj` files contain compiler flags (e.g., `UseMkl`). On non-Windows/CI, the project needs manual Fortran toolchain setup and MKL availability.
  - Running cases: the program reads a GOM project file (written by GiD preprocessor). `ReadMaterialData.for` expects the GOM header version string (Anura3D_v2025 etc.). Use small, shipped example GOM files when adding tests.

- Where to look for examples of key patterns:
  - Entry + control flow: `src/Anura3D.for`, `src/Kernel.for`, `src/InitialiseKernel.FOR`.
  - Material registration and UMAT examples: `src/ReadMaterialData.for`, `src/Soilmodels/Readme.md`, `src/Soilmodels/UMAT_LinearElasticity.f`.
  - Solver integration: `src/Solver/Solver.FOR`, `src/Solver/mkl_dss.f90`.
  - Visual Studio/Intel Fortran flags: `src/VS/Anura3D.vfproj`, `src/VS/CreateCPS.vfproj`.

- Small actionable heuristics for edits:
  - When adding a new material model: follow `Soilmodels/Readme.md` steps exactly — implement a `module`, provide `ESM_*` interface, and add `use` lines in `ExternalSoilModel.for` and `ReadMaterialData.for`.
  - Avoid touching Intel/MKL headers or VS project files unless asked; these are fragile and platform-specific.
  - Preserve fixed-form Fortran formatting; don't reflow files unless the change is strictly local.

- Tests & validation to recommend to developers:
  - Build with Visual Studio + Intel Fortran on Windows using `src/VS/Anura3D.sln`. Run a small GOM case (use any problemtype in `src/GiD_Problemtype/Anura3D_2025.gid` or the `Input_Specification/template_v2024.CPS_001`) to validate IO and material parsing.
  - For unit-like checks: add small example UMAT that returns trivial stress update; register it and run ReadMaterialData to assert `MatParams%ESM_POINTER` is set.
  - **APIC/FLIP blend feature**: Test velocity transfer by running a rotating disk or column collapse case. Adjust `CalParams%APIFLIPBlendFactor` (0-1) in `ReadCalculationData.FOR` to control PIC/FLIP blend (default 0.99). See `APIC_FLIP_IMPLEMENTATION.md` for details.

If any part of this file is unclear or you'd like more details (build scripts, known example GOM inputs, or CI hints), tell me which area to expand and I'll update the file.
