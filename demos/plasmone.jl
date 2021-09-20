using Plots

using simulator.fresnelltools
using simulator.matrixcore
using simulator.drudemetals
using simulator.utilities


angles = 0:1e-3:π/2
θ = π / 4
νs = 2e15:5e12:10e15
ν = 5e14


n_silver(ν) = c_0 / ν * (0.134 + 7.4im) - 0.05 - 0.5im
# n_silver(ν) = 0.05 + 4im
n_fiber = 1
n_air = 1.4


function system(ν, θ)
    _, interface, θ2 = FresnellBoundrary(n_fiber, n_silver(ν), θ)
    bulk = FresnellSlab(n_silver(ν), ν/c_0, 30e-9)
    _, interface2, _ = FresnellBoundrary(n_silver(ν), n_air, θ2)
    CascadeScattering([interface, bulk, interface2])
end


R(ν, θ) = map((x) -> x>1 ? 1 : x, (abs((system(ν, θ) * [1, 0])[2])^2))


resultsmatrix = Matrix{Number}(undef, length(angles), length(νs))
for (i, θ) in enumerate(angles)
    for (j, ν) in enumerate(νs)
        resultsmatrix[i, j] = R(ν, θ)
    end
end

heatmap(angles, νs, resultsmatrix, title = "TM boundrary reflectance")
