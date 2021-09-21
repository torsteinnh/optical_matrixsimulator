using Plots
plotly()

using simulator.fresnelltools
using simulator.matrixcore
using simulator.drudemetals
using simulator.utilities
using simulator.materials


angles = 0:1e-3:π/2
θ = π / 4
λs = 300e-9:1e-9:800e-9
λ = 750e-9


n_silver =  LoadMaterial("materials/Au.csv")
n_fiber = 1.4
n_air = 1
d = 20e-9


function system(λ, θ)
    _, interface, θ2 = FresnellBoundrary(n_fiber, n_silver(λ), θ)
    bulk = FresnellSlab(n_silver(λ), 2*π / λ, d / cos(θ2))
    _, interface2, _ = FresnellBoundrary(n_silver(λ), n_air, θ2)
    CascadeScattering([interface, bulk, interface2])
end


R(λ, θ) = map((x) -> x > 1 ? x : x, abs((system(λ, θ) * [1, 0])[2])^2)


Rs = R.(λ, angles)
plot(angles, Rs)
