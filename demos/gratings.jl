using Plots
import PhysicalConstants.CODATA2018: c_0
using Unitful

include("../src/simulator.jl")
using simulator.fresnelltools

n_fiber = 1.4
n_grating = 1

grating_width = 10e-9
gap_width = 400e-9

output(ν) = Grating(n_grating, n_fiber, grating_width, gap_width, 100, 2*π*ν/ustrip(c_0)) * [1, 0]


frequencies = 4e14:1e10:7e14
results = map(ν -> abs.(output(ν)).^2, frequencies)
transmission = map(pair -> pair[1], results)
reflection = map(pair -> pair[2], results)

plot(frequencies, transmission, label="transmission")
plot!(frequencies, reflection, label="reflection")
