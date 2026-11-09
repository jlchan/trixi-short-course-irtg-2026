using OrdinaryDiffEqSSPRK
using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

# Three-mode ("treble-mode") Richtmyer-Meshkov instability: a Mach 2 shock in air hits a
# sinusoidally perturbed air/helium interface. The setup follows Section 3 of
# - Pei-Guang Zhang, Jian-Ping Wang (2013)
#   A newly improved WENO scheme and its application to the simulation of
#   Richtmyer-Meshkov instability.
#   [DOI: 10.1016/j.proeng.2013.08.023](https://doi.org/10.1016/j.proeng.2013.08.023)
# The variables are nondimensionalized by the density and pressure of air at rest and by
# the velocity scale `sqrt(p0 / rho0)`. The flow is modeled as a single component: helium
# is only distinguished from air by its lower density, and both use the adiabatic index 1.4.
#
# The post-shock state is the exact Rankine-Hugoniot state of a Mach 2 shock, which the
# paper rounds to (2.67, 1.48, 0, 4.5). The interface is located at
# `x = 0.1 + 0.008 * cos(60 * pi * y)`, which matches the initial data shown in Fig. 3 of
# the paper.
function initial_condition_richtmyer_meshkov(x, t,
                                             equations::CompressibleEulerEquations2D)
    x_shock = 0.06
    x_interface = 0.1 + 0.008 * cos(60 * pi * x[2])

    if x[1] < x_shock # post-shock air
        rho, v1, v2, p = (8 / 3, 5 / 4 * sqrt(equations.gamma), 0.0, 9 / 2)
    elseif x[1] < x_interface # pre-shock air
        rho, v1, v2, p = (1.0, 0.0, 0.0, 1.0)
    else # helium
        rho, v1, v2, p = (0.138, 0.0, 0.0, 1.0)
    end

    return prim2cons(SVector(rho, v1, v2, p), equations)
end
initial_condition = initial_condition_richtmyer_meshkov

# Entropy stable DGSEM: the entropy conservative volume flux `flux_ranocha` is used with
# flux differencing, and the dissipative surface flux makes the scheme entropy stable.
surface_flux = flux_lax_friedrichs
volume_flux = flux_ranocha
polydeg = 3
basis = LobattoLegendreBasis(polydeg)

# Shock capturing: in troubled elements, the DG volume integral is blended with a
# first-order finite volume scheme on the subcells defined by the LGL nodes.
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 0.5,
                                         alpha_min = 0.001,
                                         alpha_smooth = true,
                                         variable = density_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

# A coarse mesh so that the simulation runs in a few seconds. The paper uses 1200 x 200
# grid points, i.e., `trees_per_dimension = (300, 50)` with `polydeg = 3`.
trees_per_dimension = (72, 12)
coordinates_min = (0.0, 0.0)
coordinates_max = (0.6, 0.1)
mesh = P4estMesh(trees_per_dimension; polydeg, coordinates_min, coordinates_max,
                 periodicity = (false, true))

# The transmitted shock leaves the domain through the right boundary at t ≈ 0.12.
# Imposing the initial (pre-shock) helium state there lets it leave the domain, since the
# boundary Riemann problem then consists of the outgoing shock.
boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = (; x_neg = boundary_condition,
                       x_pos = boundary_condition)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

# The final time t = 0.6 ms of the paper corresponds to t ≈ 0.174 for lengths in meters
# and the velocity scale `sqrt(R * T0) ≈ 290 m/s` of air at T0 = 293 K.
tspan = (0.0, 0.18)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43();
            adaptive = false, dt = 1, # overwritten by the `stepsize_callback`
            callback = callbacks);

###############################################################################
# visualize the density at the final time

using Plots
pd = PlotData2D(sol)
plot(pd["rho"], aspect_ratio = 1, size = (900, 250),
     title = "Density at t = $(round(sol.t[end], digits = 3))")
