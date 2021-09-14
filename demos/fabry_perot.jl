using Plots
include("../src/utilities.jl")
using .utilities

include("../src/simulator.jl")
using simulator.matrixcore
using simulator.fresnelltools

mirror1_r = 0.6
mirror2_r = 0.6
distance = 3e-4

mirror1 = [
    1 - mirror1_r mirror1_r;
    mirror1_r     1 - mirror1_r
]
mirror2 = [
    1 - mirror2_r mirror2_r;
    mirror2_r     1 - mirror2_r
]

gap(ν) = FresnellSlab(1, 2*π*ν/c_0, distance)

transmittance(ν) = CascadeScattering([mirror1, gap(ν), mirror2])[1, 1]
golden(ν) = (1 - mirror1_r) * (1 - mirror2_r) * ℯ^(-im * distance * 2*π*ν / c_0) / (1 - mirror1_r * mirror2_r * ℯ^(-2im * distance * 2*π*ν / c_0))

xs = 1e9:1e8:2e12
transmission = map(ν -> abs(transmittance(ν))^2, xs)
transmission_g = map(ν -> abs(golden(ν))^2, xs)

reflection = 1 .- transmission
reflection_g = 1 .- transmission_g

plot(xs, transmission, label="transmission")
plot!(xs, reflection, label="reflection")
plot!(xs, transmission_g, label="transmission_g")
plot!(xs, reflection_g, label="reflection_g")
