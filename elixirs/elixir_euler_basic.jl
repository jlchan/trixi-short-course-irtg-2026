# The anatomy of a Trixi.jl elixir: a minimal entropy stable DGSEM setup for the 2D
# compressible Euler equations. Visualize the result with `using Plots; plot(sol)`.
using OrdinaryDiffEqLowStorageRK
using Trixi

equations = CompressibleEulerEquations2D(1.4)
initial_condition = initial_condition_weak_blast_wave

solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs,
               volume_integral = VolumeIntegralFluxDifferencing(flux_ranocha))

mesh = TreeMesh((-2.0, -2.0), (2.0, 2.0),
                initial_refinement_level = 5, periodicity = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition_periodic)
ode = semidiscretize(semi, (0.0, 0.4))

callbacks = CallbackSet(SummaryCallback(),
                        AnalysisCallback(semi, interval = 100),
                        StepsizeCallback(cfl = 1.0))

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # overwritten by the `StepsizeCallback`
            callback = callbacks, ode_default_options()...)
