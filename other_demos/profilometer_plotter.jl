using Plots
using DelimitedFiles

save = false


if save
    gr()
else
    plotly()
end


function LoadStylus(path, flattening)
    data = readdlm(path, ',', skipstart=34)[:, 1:2]

    correction = data[flattening, 2]
    for i in 1:length(data[:, 2])
        data[i, 2] = data[i, 2] - correction * i / flattening
    end

    data[:, 1], data[:, 2]
end


stylusfig = plot(title="Stylus measurements", xlabel="distance [μm]", ylabel="height [nm]", ticks=:native)
for sample in ["5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]
    xs, ys = LoadStylus("../NanoLab_process/20220405/S" * sample * ".csv", 16000)

    plot!(stylusfig, xs, ys, label=sample)
end
display(stylusfig)
if save savefig(stylusfig, "../../rapport/figures/profiles/stylus.pdf") end


stylusfig_contested = plot(title="Stylus measurements", xlabel="distance [μm]", ylabel="height [nm]", ticks=:native)
for sample in [
        "20220405/S5", "20220424/5A",
        "20220405/S8", "20220424/8A",
        "20220405/S9", "20220424/9A",
        "20220405/S10", "20220424/10A",
        "20220405/S12", "20220424/12A",
        "20220405/S13", "20220424/13A",
        "20220405/S14", "20220424/14A",
        "20220405/S18", "20220424/18A",
        "20220405/S19", "20220424/19A",
        "20220405/S20", "20220424/20A",
        ]
    xs, ys = LoadStylus("../NanoLab_process/" * sample * ".csv", 15500)

    plot!(stylusfig_contested, xs, ys, label=sample)
end
display(stylusfig_contested)
if save savefig(stylusfig_contested, "../../rapport/figures/profiles/stylus.pdf") end


stylusfig_extra = plot(title="Stylus measurements", xlabel="distance [μm]", ylabel="height [nm]", ticks=:native)
for sample in ["21", "22", "23", "24"]
    xs, ys = LoadStylus("../NanoLab_process/20220505/" * sample * "A.csv", 8000)

    plot!(stylusfig_extra, xs, ys, label=sample)
end
display(stylusfig_extra)
if save savefig(stylusfig_extra, "../../rapport/figures/profiles/stylus.pdf") end