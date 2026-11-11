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

# A coarse base mesh, which is adaptively refined near the shocks and the interface (see the
# `AMRCallback` below). The paper uses 1200 x 200 grid points, i.e., about
# `trees_per_dimension = (300, 50)` with `polydeg = 3`. Two levels of refinement of the base
# mesh give a comparable resolution of 288 x 48 elements.
trees_per_dimension = (72, 12)
coordinates_min = (0.0, 0.0)
coordinates_max = (0.6, 0.1)
mesh = P4estMesh(trees_per_dimension; polydeg, coordinates_min, coordinates_max,
                 periodicity = (false, true))

# Custom implementation of the characteristic upwind boundary condition 
# for the 2D compressible Euler equations
struct BoundaryConditionCharacteristicUpwind{B}
    boundary_value_function::B
end

# `normal_direction` is the outward normal scaled by the surface Jacobian
@inline function (bc::BoundaryConditionCharacteristicUpwind)(u,
                                                             normal_direction::AbstractVector,
                                                             x, t, surface_flux,
                                                             equations::CompressibleEulerEquations2D)
    norm_ = sqrt(normal_direction[1]^2 + normal_direction[2]^2)
    n1, n2 = normal_direction / norm_
    rho, v1, v2, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)
    v_n = v1 * n1 + v2 * n2
    H = (u[4] + p) / rho
    # eigenvalues and right eigenvectors of the flux Jacobian in the unit normal
    # direction for the (acoustic -, shear, entropy, acoustic +) waves
    lambda = SVector(v_n - c, v_n, v_n, v_n + c)
    R = hcat(SVector(1, v1 - c * n1, v2 - c * n2, H - c * v_n),
             SVector(0, -n2, n1, v2 * n1 - v1 * n2),
             SVector(1, v1, v2, 0.5f0 * (v1^2 + v2^2)),
             SVector(1, v1 + c * n1, v2 + c * n2, H + c * v_n))
    u_ext = bc.boundary_value_function(x, t, equations)
    return flux(u, normal_direction, equations) +
           norm_ * (R * (min.(lambda, 0) .* (R \ (u_ext - u))))
end

# The post-shock air state is imposed at the left (inflow) boundary. The transmitted shock
# leaves the domain through the right boundary at t ≈ 0.12. There, the characteristic upwind
# boundary condition lets outgoing waves leave the domain and only takes incoming waves from
# the initial (pre-shock) helium state.
boundary_conditions = (; x_neg = BoundaryConditionDirichlet(initial_condition),
                         x_pos = BoundaryConditionCharacteristicUpwind(initial_condition))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

# The final time t = 0.6 ms of the paper corresponds to t ≈ 0.174 for lengths in meters
# and the velocity scale `sqrt(R * T0) ≈ 290 m/s` of air at T0 = 293 K.
tspan = (0.0, 0.18)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

# Adaptive mesh refinement: the Löhner indicator uses a normalized second difference 
# of the density). Elements are refined up to `max_level` above the base mesh and 
# coarsened again once the features have passed.
amr_indicator = IndicatorLöhner(semi, variable = Trixi.density)
amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level = 0,
                                      med_level = 1, med_threshold = 0.05,
                                      max_level = 2, max_threshold = 0.1)
amr_callback = AMRCallback(semi, amr_controller,
                           interval = 5,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true)

stepsize_callback = StepsizeCallback(cfl = 1.0)

# The `stepsize_callback` comes after the `amr_callback` so that the time step is computed
# on the adapted mesh.
callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        amr_callback,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43();
            adaptive = false, dt = 1, # overwritten by the `stepsize_callback`
            callback = callbacks, 
            ode_default_options()...);

###############################################################################
# visualize the density and the adapted mesh at the final time

using Plots
pd = PlotData2D(sol)
p_rho = plot(pd["rho"], aspect_ratio = 1,
             title = "Density at t = $(round(sol.t[end], digits = 3))")
p_mesh = plot(getmesh(pd), aspect_ratio = 1, linewidth = 0.5,
              xlims = (0.0, 0.6), ylims = (0.0, 0.1), title = "Adapted mesh")
p = plot(p_rho, p_mesh, layout = (2, 1), size = (900, 500))
savefig(p, "sol.png")