using Plots
plotly()

using simulator.fresnelltools
using simulator.materials.spesifics
using simulator.matrixcore


θs = 0:1e-3:π/2
θ = 0.806
λs = 300e-9:1e-9:1200e-9
λ = 800e-9
ds = 1e-9:1e-9:100e-9
d = 53e-9

n_palladium = Pd
n_fiber = 1.4
n_air = 1


function structure(λ, θ, d_palladium)
    interface_s_pd_te, interface_s_pd_tm, θpd = FresnellBoundrary(n_fiber, n_palladium(λ), θ)
    bulk_pd = FresnellSlab(n_palladium(λ), 2*π/λ, d_palladium, θpd)
    interface_pd_a_te, interface_pd_a_tm, θa = FresnellBoundrary(n_palladium(λ), n_air, θpd)

    te_system = CascadeScattering([interface_s_pd_te, bulk_pd, interface_pd_a_te])
    tm_system = CascadeScattering([interface_s_pd_tm, bulk_pd, interface_pd_a_tm])

    te_system, tm_system
end

R_te(λ, θ, d) = abs((structure(λ, θ, d)[1] * [1, 0])[2])^2
R_tm(λ, θ, d) = abs((structure(λ, θ, d)[2] * [1, 0])[2])^2
R_ratio(λ, θ, d) = R_te(λ, θ, d) / R_tm(λ, θ, d)


# plot(λs, R_ratio.(λs, θ, d), label="θ=0.856", title="TE Reflection / TM Reflection at d=53nm")
# plot(θs, (θ) -> (R_ratio.(λ, θ + 1e-9im, d)))
# plot(ds, R_ratio.(λ, θ, ds))

plot(λs, R_te.(λs, θ, d), label="te", title="Reflection at silica-palladium-air")
plot!(λs, R_tm.(λs, θ, d), label="tm", xlabel="incident angle", ylabel="reflected power")