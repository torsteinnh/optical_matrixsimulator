using Plots
plotly()

using simulator.fresnelltools
using simulator.materials.spesifics
using simulator.matrixcore


θs = 0:1e-3:π/2
θ = 0.806
λs = 350e-9:1e-9:1800e-9
λ = 800e-9

d_ags = 10e-9:1e-9:200e-9
d_ag = 50e-9
d_sis = 10e-9:1e-9:200e-9
d_si = 70e-9
d_pds = 10e-9:1e-9:200e-9
d_pd = 5e-9

n_fiber = SiO2_core_Sellmeier
n_ag = Ag
n_silicafilm = SiO2_thinfilm_Ciprian
n_palladium(λ) = Pd(λ) - 0im
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

R_te(λ, θ, d_ag, d_si, d_pd) = abs((structure(λ, θ, d_ag, d_si, d_pd)[1] * [1, 0])[2])^2
R_tm(λ, θ, d_ag, d_si, d_pd) = abs((structure(λ, θ, d_ag, d_si, d_pd)[2] * [1, 0])[2])^2
R_ratio(λ, θ, d_ag, d_si, d_pd) = R_te(λ, θ, d_ag, d_si, d_pd) / R_tm(λ, θ, d_ag, d_si, d_pd)


# plot(λs, R_ratio.(λs, θ, d_ag, d_si, d_pd), label="θ=0.856", title="TE Reflection / TM Reflection at d=53nm")
# plot(θs, (θ) -> (R_ratio.(λ, θ + 1e-9im, d_ag, d_si, d_pd)))
# plot(d_ags, R_ratio.(λ, θ, d_ags, d_si, d_pd))
# plot(d_sis, R_ratio.(λ, θ, d_ag, d_sis, d_pd))
# plot(d_pds, R_ratio.(λ, θ, d_ag, d_si, d_pds))

plot(λs, R_te.(λs, θ, d_ag, d_si, d_pd), label="te", title="Reflection at silica-palladium-air")
plot!(λs, R_tm.(λs, θ, d_ag, d_si, d_pd), label="tm", xlabel="incident angle", ylabel="reflected power")
