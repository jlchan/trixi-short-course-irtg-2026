# Trixi.jl short course (IRTG 2026)

Material for a short course on structure-preserving high-order discontinuous Galerkin
methods with [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/jlchan/trixi-short-course-irtg-2026?quickstart=1)

## Getting started in GitHub Codespaces

1. Click the button above or, on the GitHub page of this repository, click
   **Code → Codespaces → Create codespace on main**.
2. Wait until the setup has finished. It installs Julia and precompiles all packages,
   which can take several minutes for a new codespace.
3. Start a Julia REPL (command palette: **Julia: Start REPL**) and run
   ```julia
   include("elixirs/elixir_euler_richtmyer_meshkov.jl")
   ```
   The first run takes longer since Julia compiles the code; later runs are much faster.

Codespaces count against the monthly usage quota of your GitHub account. Stop or delete
your codespace when you no longer need it.

## Running locally

1. Install Julia 1.12, e.g., via [juliaup](https://github.com/JuliaLang/juliaup).
2. Clone this repository and install the packages:
   ```shell
   git clone https://github.com/jlchan/trixi-short-course-irtg-2026.git
   cd trixi-short-course-irtg-2026
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
3. Start Julia with `julia --project=.` and run the elixir as above.

## Contents

- [`elixirs/elixir_advection_basic.jl`](elixirs/elixir_advection_basic.jl):
  The simplest Trixi.jl setup, 2D linear advection with periodic boundaries, which shows
  the same anatomy as `elixir_euler_basic.jl` with a single DGSEM solver line.
- [`elixirs/elixir_euler_basic.jl`](elixirs/elixir_euler_basic.jl):
  A minimal elixir that shows the anatomy of a Trixi.jl simulation (equations, solver,
  mesh, semidiscretization, callbacks, time integration) for a weak blast wave of the
  compressible Euler equations, using entropy stable flux differencing on a `TreeMesh`.
- [`elixirs/elixir_euler_richtmyer_meshkov.jl`](elixirs/elixir_euler_richtmyer_meshkov.jl):
  Richtmyer-Meshkov instability for the compressible Euler equations, using entropy
  stable flux differencing with shock capturing on a `P4estMesh`.

See [`elixirs/README.md`](elixirs/README.md) for things to try with the Richtmyer-Meshkov
elixir.
