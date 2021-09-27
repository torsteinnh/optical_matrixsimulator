using Plots
plotly()

using simulator.fresnelltools
using simulator.materials


θs = 0:1e-3:π/2
θ = 0.856
λs = 300e-9:1e-9:1200e-9
λ = 797e-9
ds = 1e-9:1e-9:100e-9
d = 53e-9


n_gold = LoadMaterial("materials/Au.csv")
n_fiber = 1.4
n_air = 1


R_te(λ, θ, d) = abs((ThreeLayerSystem(n_fiber, n_gold, n_air, λ, θ, d)[1] * [1, 0])[2])^2
R_tm(λ, θ, d) = abs((ThreeLayerSystem(n_fiber, n_gold, n_air, λ, θ, d)[2] * [1, 0])[2])^2
R_ratio(λ, θ, d) = R_te(λ, θ, d) / R_tm(λ, θ, d)


# plot(λs, R_ratio.(λs, θ, d), label="θ=0.856", title="TE Reflection / TM Reflection at d=53nm")
# plot(θs, R_ratio.(λ, θs, d))
# plot(ds, R_ratio.(λ, θ, ds))

# plot(θs, R_te.(λ, θs, d), label="te", title="Reflection at λ=797nm, d=90nm")
# plot!(θs, R_tm.(λ, θs, d), label="tm")
