using Plots
plotly()

using simulator.fresnelltools
using simulator.materials.spesifics


θs = 0:1e-3:π/2
θ = 0.82
λs = 300e-9:1e-9:1200e-9
λ = 797e-9
ds = 1e-9:1e-9:100e-9
d = 53e-9


n_gold = Au
n_fiber = 1.4
n_air = 1


R_te(λ, θ, d) = abs((ThreeLayerSystem(n_fiber, n_gold, n_air, λ, θ, d)[1] * [1, 0])[2])^2
R_tm(λ, θ, d) = abs((ThreeLayerSystem(n_fiber, n_gold, n_air, λ, θ, d)[2] * [1, 0])[2])^2
R_ratio(λ, θ, d) = R_te(λ, θ, d) / R_tm(λ, θ, d)


function Inside(λ, θ, d)
    border_1_te, border_1_tm, θg = FresnellBoundrary(n_fiber, n_gold(λ), θ)
    border_2_te, border_2_tm, θa = FresnellBoundrary(n_gold(λ), n_air, θg)

    bulk_fiber(d) = FresnellSlab(n_fiber, 2*π/λ, d, θ)
    bulk_gold(d) = FresnellSlab(n_gold(λ), 2*π/λ, d, θg)
    bulk_air(d) = FresnellSlab(n_air, 2*π/λ, d, θa)

    border_1_te_f(d) = border_1_te
    border_1_tm_f(d) = border_1_tm
    border_2_te_f(d) = border_2_te
    border_2_tm_f(d) = border_2_tm

    search = (Up, Un) -> abs(Up + Un)^2

    u_te, d_te = Interrogator(
        [bulk_fiber, border_1_te_f, bulk_gold, border_2_te_f, bulk_air],
        [20e-9, 0, d, 0, 20e-9],
        1e-10,
        search)

    u_tm, d_tm = Interrogator(
        [bulk_fiber, border_1_tm_f, bulk_gold, border_2_tm_f, bulk_air],
        [20e-9, 0, d, 0, 20e-9],
        1e-10,
        search)

    d_te, u_te, u_tm
end


# plot(λs, R_ratio.(λs, θ, d), label="θ=0.856", title="TE Reflection / TM Reflection at d=53nm")
# plot(θs, (θ) -> (R_ratio.(λ, θ + 1e-9im, d)))
# plot(ds, R_ratio.(λ, θ, ds))

# plot(λs, R_te.(λs, θ, d), label="te", title="Reflection at θ=0.82, d=53nm", legend=:bottomright)
# plot!(λs, R_tm.(λs, θ, d), label="tm", xlabel="wavelength", ylabel="reflected power")


ds, u_te, u_tm = Inside(650e-9, θ, d)
ds = ds .* 1e9
_, u_tep, u_tmp = Inside(786e-9, θ, d)

plot(ds, u_te, label="te normal", title="Gold thin film, power distribution",
    xlabel="distance in nm", ylabel="all propagating power",
    legend=:topleft)
plot!(ds, u_tm, label="tm normal")
plot!(ds, u_tep, label="te plasmone")
plot!(ds, u_tmp, label="tm plasmone")
