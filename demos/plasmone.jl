using Plots
plotly()

using simulator.fresnelltools
using simulator.matrixcore
using simulator.drudemetals
using simulator.utilities
using simulator.materials


θs = 0:1e-3:π/2
θ = 0.7
λs = 300e-9:1e-9:1200e-9
λ = 797e-9


n_silver =  LoadMaterial("materials/Au.csv")
n_fiber = 1.4
n_air = 1
d = 20e-9


function system(λ, θ)
    interface1_te, interface1_tm, θ2 = FresnellBoundrary(n_fiber, n_silver(λ), θ)
    bulk = FresnellSlab(n_silver(λ), 2*π / λ, d, θ2)
    interface2_te, interface2_tm, _ = FresnellBoundrary(n_silver(λ), n_air, θ2)
    te_system = CascadeScattering([interface1_te, bulk, interface2_te])
    tm_system = CascadeScattering([interface1_tm, bulk, interface2_tm])

    te_system, tm_system
end


R_te(λ, θ) = abs((system(λ, θ)[1] * [1, 0])[2])^2
R_tm(λ, θ) = abs((system(λ, θ)[2] * [1, 0])[2])^2
R_ratio(λ, θ) = R_te(λ, θ) / R_tm(λ, θ)


# plot(λs, R_ratio.(λs, θ))
plot(θs, R_ratio.(λ, θs))
