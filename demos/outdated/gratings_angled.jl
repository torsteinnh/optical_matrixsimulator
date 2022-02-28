using Plots

using simulator.fresnelltools
include("../../src/utilities.jl")
using .utilities


n_normal = 1.4
n_spechial = 1

grating_width = 10e-9
gap_width = 400e-9

frequencies = 4e14:1e11:7e14
angles = 0:1e-2:(π/2)

layers = 100


gratings(ν, θ) = Grating(n_spechial, n_normal, grating_width, gap_width, layers, 2*π*ν/c_0, θ)
get_r_te(ν, θ) = abs((gratings(ν, θ)[1] * [1, 0])[2])^2
get_r_tm(ν, θ) = abs((gratings(ν, θ)[2] * [1, 0])[2])^2


function get_data(angles, frequencies)
    te_data = Matrix{Real}(undef, length(angles), length(frequencies))
    tm_data = Matrix{Real}(undef, length(angles), length(frequencies))
    for (i, θ) in enumerate(angles)
        for (j, ν) in enumerate(frequencies)
            te_data[i, j] = get_r_te(ν, θ)
            tm_data[i, j] = get_r_tm(ν, θ)
        end
    end
    te_data, tm_data
end

(te_data, tm_data) = get_data(angles, frequencies)
heatmap(angles, frequencies, te_data, xlabel="angles", ylabel="frequencies", title="TE modes")
heatmap(angles, frequencies, tm_data, xlabel="angles", ylabel="frequencies", title="TM modes")
