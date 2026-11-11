# Things to try with the Richtmyer-Meshkov elixir

`trixi_include` runs an elixir with some of its assignments or keyword arguments replaced,
without editing the file. In a Julia REPL started with `julia --project=.` in the repository
root (or in the codespace):

```julia
using Trixi
elixir = "elixirs/elixir_euler_richtmyer_meshkov.jl"

# a uniform mesh at the resolution of the finest AMR level, without AMR
trixi_include(elixir; trees_per_dimension = (288, 48), amr_callback = TrivialCallback())

# weaker shock capturing (the elixir uses alpha_max = 0.5); what happens with alpha_max = 0?
trixi_include(elixir; alpha_max = 0.1)

# error-based adaptive time steps of OrdinaryDiffEq.jl instead of the CFL condition
trixi_include(elixir; stepsize_callback = TrivialCallback(), adaptive = true, dt = 1.0e-4)

# reshock: a wall at the right boundary reflects the transmitted shock into the mixing zone
trixi_include(elixir; x_pos = boundary_condition_slip_wall, tspan = (0.0, 0.3))

# a higher polynomial degree, which needs a smaller CFL number
trixi_include(elixir; polydeg = 5, cfl = 0.5)

# another surface flux and one more level of refinement
trixi_include(elixir; surface_flux = flux_hllc, max_level = 3)
```

The following changes need edits of the elixir (work on a copy, e.g.,
`cp elixirs/elixir_euler_richtmyer_meshkov.jl my_elixir.jl`).

**Zhang-Shu positivity limiting.** Define the limiter after `semi`, pass
`limiter! = positivity_limiter` to the `AMRCallback`, and use
`SSPRK43(stage_limiter! = positivity_limiter)` in `solve`:
```julia
positivity_limiter = PositivityPreservingLimiterZhangShu(thresholds = (1.0e-6, 1.0e-6),
                                                         variables = (Trixi.density, pressure))
```

**Liu-Zhang positivity limiting.** In addition to the Zhang-Shu limiter above (which stays in
the `AMRCallback`), define the global limiter and use
`SSPRK43(stage_limiter! = global_limiter, step_limiter! = global_limiter)`:
```julia
global_limiter = PositivityPreservingLimiterLiuZhang(positivity_limiter, semi)
```

**Entropy correction instead of shock capturing**, together with the Liu-Zhang limiter:
```julia
volume_integral = VolumeIntegralAdaptive(IndicatorEntropyCorrection(equations, basis;
                                                                    scaling = 2),
                                         VolumeIntegralWeakForm(),
                                         VolumeIntegralPureLGLFiniteVolume(surface_flux))
```

**Subcell limiting instead of shock capturing**, with `StepsizeCallback(cfl = 0.5)` and
Trixi.jl's own SSP time integrator:
```julia
limiter_idp = SubcellLimiterIDP(equations, basis;
                                positivity_variables_cons = ["rho"],
                                positivity_variables_nonlinear = [pressure])
volume_integral = VolumeIntegralSubcellLimiting(limiter_idp;
                                                volume_flux_dg = volume_flux,
                                                volume_flux_fv = surface_flux)
# ...
sol = Trixi.solve(ode, Trixi.SimpleSSPRK33(stage_callbacks = (SubcellLimiterIDPCorrection(),));
                  dt = 1.0, callback = callbacks, ode_default_options()...);
```

**A heavy gas instead of helium.** In `initial_condition_richtmyer_meshkov`, replace the
helium density `0.138` by `5.0` (about the density ratio of SF6 and air).

<details>
<summary>What we observed (Trixi.jl v0.17.21)</summary>

- The uniform 288 x 48 mesh takes about 2.5 times as long as the AMR run.
- `alpha_max = 0.1` and `0.02` run; with `alpha_max = 0` (no shock capturing), the run fails at
  t ≈ 0.016: entropy stability alone does not keep density and pressure positive.
- With the default shock capturing, the Zhang-Shu and Liu-Zhang limiters are never active.
  Without shock capturing, they keep the run going, but the density oscillates and the time
  step becomes very small.
- Entropy correction without a limiter fails at t ≈ 0.014. With the Liu-Zhang limiter,
  `scaling = 2` runs (about 1.8 times the runtime of the original elixir); `scaling = 1` also runs,
  but about 10 times longer and with strong under- and overshoots.
- Subcell limiting with positivity bounds needs `cfl = 0.5`; with `cfl = 1`, the low-order
  method itself loses positivity.
- Error-based time stepping takes about 15% more steps than the CFL-based step size here, and
  rejects about one step in six.
- `polydeg = 5` fails with `cfl = 1.0`, with or without a positivity limiter, and runs with
  `cfl = 0.5`.
- With the heavy gas, the interface does not invert its phase (with helium it does).

</details>
