using Plots
plotly()

using simulator.fresnelltools
using simulator.materials.spesifics
using simulator.matrixcore


θs = 0:1e-3:π/2
θ = 0.806
λs = 350e-9:1e-9:1800e-9
λ = 1434e-9

d_ags = 10e-9:1e-9:200e-9
d_ag = 50e-9
d_sis = 10e-9:1e-9:200e-9
d_si = 20e-9
d_pds = 10e-9:1e-9:200e-9
d_pd = 5e-9

n_fiber = SiO2_core_Sellmeier
n_ag = Ag
n_silicafilm = SiO2_thinfilm_Ciprian
n_palladium(λ) = Pd_Werner(λ) - 0im
n_air = 1


function structure(λ, θ, d_ag, d_si, d_pd)
    interface_fiber_ag_te, interface_fiber_ag_tm, θag = FresnellBoundrary(n_fiber(λ), n_ag(λ), θ)
    bulk_ag = FresnellSlab(n_ag(λ), 2*π/λ, d_ag, θag)
    interface_ag_si_te, interface_ag_si_tm, θsi = FresnellBoundrary(n_ag(λ), n_silicafilm(λ), θag)
    bulk_si = FresnellSlab(n_silicafilm(λ), 2*π/λ, d_si, θsi)
    interface_si_pd_te, interface_si_pd_tm, θpd = FresnellBoundrary(n_silicafilm(λ), n_palladium(λ), θsi)
    bulk_pd = FresnellSlab(n_palladium(λ), 2*π/λ, d_pd, θpd)
    interface_pd_air_te, interface_pd_air_tm, _ = FresnellBoundrary(n_palladium(λ), n_air, θpd)

    te_system = CascadeScattering([interface_fiber_ag_te, bulk_ag, interface_ag_si_te, bulk_si, interface_si_pd_te, bulk_pd, interface_pd_air_te])
    tm_system = CascadeScattering([interface_fiber_ag_tm, bulk_ag, interface_ag_si_tm, bulk_si, interface_si_pd_tm, bulk_pd, interface_pd_air_tm])

    te_system, tm_system
end

function structure_inside(λ, θ, d_ag, d_si, d_pd)
    bulk_fiber(d) = FresnellSlab(n_fiber(λ), 2*π/λ, d, θ)
    interface_fiber_ag_te, interface_fiber_ag_tm, θ_ag = FresnellBoundrary(n_fiber(λ), n_ag(λ), θ)
    interface_fiber_ag_tm_f(d) = interface_fiber_ag_tm
    interface_fiber_ag_te_f(d) = interface_fiber_ag_te
    bulk_ag(d) = FresnellSlab(n_ag(λ), 2*π/λ, d, θ_ag)
    interface_ag_si_te, interface_ag_si_tm, θ_si = FresnellBoundrary(n_ag(λ), n_silicafilm(λ), θ_ag)
    interface_ag_si_tm_f(d) = interface_ag_si_tm
    interface_ag_si_te_f(d) = interface_ag_si_te
    bulk_si(d) = FresnellSlab(n_silicafilm(λ), 2*π/λ, d, θ_si)

    # interface_si_pd_te, interface_si_pd_tm, θ_pd = FresnellBoundrary(n_silicafilm(λ), n_palladium(λ), θ_si)
    # interface_si_pd_te_f(d) = interface_si_pd_te
    # interface_si_pd_tm_f(d) = interface_si_pd_tm
    # bulk_pd(d) = FresnellSlab(n_palladium(λ), 2*π/λ, d, θ_pd)
    # interface_pd_air_te, interface_pd_air_tm, θ_air = FresnellBoundrary(n_palladium(λ), n_air, θ_air)
    # interface_pd_air_te_f(d) = interface_pd_air_te
    # interface_pd_air_tm_f(d) = interface_pd_air_tm
    # bulk_air(d) = FresnellSlab(n_air, 2*π/λ, d, θ_air)
    
    interface_si_air_te, interface_si_air_tm, θ_air = FresnellBoundrary(n_silicafilm(λ), n_air, θ_si)
    interface_si_air_te_f(d) = interface_si_air_te
    interface_si_air_tm_f(d) = interface_si_air_tm
    bulk_air(d) = FresnellSlab(n_air, 2*π/λ, d, θ_air)

    interrogated = (Up, Un) -> abs(Up + Un)^2

    u_values_tm, d_values_tm = Interrogator(
        [bulk_fiber, interface_fiber_ag_tm_f, bulk_ag, interface_ag_si_tm_f, bulk_si, interface_si_air_tm_f, bulk_air],
        [10e-9, 0, d_ag, 0, d_si, 0, 10e-9], 1e-10,
        interrogated)
    u_values_te, d_values_te = Interrogator(
        [bulk_fiber, interface_fiber_ag_te_f, bulk_ag, interface_ag_si_te_f, bulk_si, interface_si_air_te_f, bulk_air],
        [10e-9, 0, d_ag, 0, d_si, 0, 10e-9], 1e-10,
        interrogated)

    u_values_te, d_values_te, u_values_tm, d_values_tm
end

R_te(λ, θ, d_ag, d_si, d_pd) = abs((structure(λ, θ, d_ag, d_si, d_pd)[1] * [1, 0])[2])^2
R_tm(λ, θ, d_ag, d_si, d_pd) = abs((structure(λ, θ, d_ag, d_si, d_pd)[2] * [1, 0])[2])^2
R_ratio(λ, θ, d_ag, d_si, d_pd) = R_te(λ, θ, d_ag, d_si, d_pd) / R_tm(λ, θ, d_ag, d_si, d_pd)

interrogated_u_te, interrogated_d_te, interrogated_u_tm, interrogated_d_tm = structure_inside(λ, θ, d_ag, d_si, d_pd)


# plot(λs, R_ratio.(λs, θ, d_ag, d_si, d_pd), label="θ=0.856", title="TE Reflection / TM Reflection at d=53nm")
# plot(θs, (θ) -> (R_ratio.(λ, θ + 1e-9im, d_ag, d_si, d_pd)))
# plot(d_ags, R_ratio.(λ, θ, d_ags, d_si, d_pd))
# plot(d_sis, R_ratio.(λ, θ, d_ag, d_sis, d_pd))
# plot(d_pds, R_ratio.(λ, θ, d_ag, d_si, d_pds))

# plot(λs, R_te.(λs, θ, d_ag, d_si, d_pd), label="te, Pd-structure, Different Pd", title="Reflection at silica-palladium-air", legend=:outertopright)
# plot!(λs, R_tm.(λs, θ, d_ag, d_si, d_pd), label="tm, Pd-structure, Different Pd", xlabel="wavelength", ylabel="reflected power")

# plot(interrogated_d_te .* 1e9, interrogated_u_te, label="te, λ=$λ", title="Sum of fields at plasmone conditions", legend=:topright, xlabel="distance in structure", ylabel="total power")
# plot!(interrogated_d_tm .* 1e9, interrogated_u_tm, label="tm λ=$λ")
