# The anatomy of a Trixi.jl elixir: a minimal DGSEM setup for 2D linear advection with
# periodic boundaries. Visualize the result with `using Plots; plot(sol)`.
using OrdinaryDiffEqLowStorageRK
using Trixi

equations = LinearScalarAdvectionEquation2D((0.2, -0.7))
initial_condition = initial_condition_convergence_test

solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

mesh = TreeMesh((-1.0, -1.0), (1.0, 1.0),
                initial_refinement_level = 4, periodicity = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition_periodic)
ode = semidiscretize(semi, (0.0, 1.0))

callbacks = CallbackSet(SummaryCallback(),
                        AnalysisCallback(semi, interval = 100),
                        StepsizeCallback(cfl = 1.6))

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # overwritten by the `StepsizeCallback`
            callback = callbacks, ode_default_options()...)
