using Plots

using simulator.fresnelltools
using simulator.matrixcore
using simulator.materials.spesifics

save = false

n_1 = SiO2_core_Sellmeier
n_2 = Au
n_3 = Pd
n_4 = l -> 1

d2 = 10e-9
d2s = 0e-9:1e-9:40e-9
d3s = 0e-9:1e-9:60e-9

λs = 900e-9:1e-9:1250e-9
λs_broadband = 400e-9:1e-9:1800e-9
θd = 45
θds = 5:1:80

description = "Palladium on $(round(d2 * 1e9, digits=2)) nm gold"


fig_path = "../rapport/figures/thickness_example/"
fig_token = replace(description, " " => "_")

if save
    gr()
else
    plotly()
end
θ1 = θd * π / 180
θs = θds .* (π / 180)

function system_single(λ::Float64, θ1::Float64, d2::Float64, d3::Float64)
    i_1_2_te, i_1_2_tm, θ2 = FresnellBoundrary(n_1(λ), n_2(λ), θ1)
    bulk_2 = FresnellSlab(n_2(λ), 2*π/λ, d2, θ2)
    i_2_3_te, i_2_3_tm, θ3 = FresnellBoundrary(n_2(λ), n_3(λ), θ2)
    bulk_3 = FresnellSlab(n_3(λ), 2*π/λ, d3, θ3)
    i_3_4_te, i_3_4_tm, _ = FresnellBoundrary(n_3(λ), n_4(λ), θ3)

    system_te = CascadeScattering([i_1_2_te, bulk_2, i_2_3_te, bulk_3, i_3_4_te])
    system_tm = CascadeScattering([i_1_2_tm, bulk_2, i_2_3_tm, bulk_3, i_3_4_tm])

    te_r = abs(system_te[2, 1])^2
    te_t = abs(system_te[1, 1])^2
    tm_r = abs(system_tm[2, 1])^2
    tm_t = abs(system_tm[1, 1])^2

    [te_r, tm_r, te_t, tm_t]
end

function get_min(λs, θ1::Float64, d2::Float64, d3::Float64)
    mins = [2.0, 2, 2, 2]
    λmins = [0.0, 0, 0, 0]

    for λ in λs
        spesifics = system_single(λ, θ1, d2, d3)
        for (i, min, spesific) in zip(1:4, mins, spesifics)
            if spesific <= min
                mins[i] = spesific
                λmins[i] = λ
            end
        end
    end

    mins, λmins
end

println("Settup complete, starting simulation")


# get_min_results = [get_min(λs, θ1, d2, d3_spesific) for d3_spesific in d3s]
# min_r_te = [item[1][1] for item in get_min_results]
# min_r_tm = [item[1][2] for item in get_min_results]
# min_λ_te = [item[2][1] for item in get_min_results]
# min_λ_tm = [item[2][2] for item in get_min_results]

## Plot plasmon minima for plasmon as a function of d3
# fig_r = plot(title="Minimum te and tm reflection for $description at θ = $θd degrees", legend=:topright, xlabel="thickness in nm", ylabel="minimum reflection coefficient", ticks=:native)
# plot!(fig_r, d3s .* 1e9, min_r_te, label="TE reflection")
# plot!(fig_r, d3s .* 1e9, min_r_tm, label="TM reflection")
# display(fig_r)
# if save savefig(fig_r, fig_path * "min_reflection_" * fig_token * ".pdf") end

## Plot frequency minima for plasmon as a function of d3
# fig_λ = plot(title="λ for minimum te and tm reflection for $description at θ = $θd degrees", legend=:topright, xlabel="thickness in nm", ylabel="λ at minimum reflection coefficient", ticks=:native)
# plot!(fig_λ, d3s .* 1e9, min_λ_te, label="TE reflection")
# plot!(fig_λ, d3s .* 1e9, min_λ_tm, label="TM reflection")
# display(fig_λ)
# if save savefig(fig_λ, fig_path * "lambda_min_reflection_" * fig_token * ".pdf") end

## Plot plasmon minima and λ for plasmon as a function of d3 for all d2s
# fig_r_all = plot(title="Minimum tm reflection for $description at θ = $θd degrees", legend=:topright, xlabel="thickness in nm", ylabel="minimum reflection coefficient", ticks=:native)
# fig_λ_all = plot(title="λ for minimum tm reflection for $description at θ = $θd degrees", legend=:topright, xlabel="thickness in nm", ylabel="λ at minimum reflection coefficient", ticks=:native)
# for d2_spesific in d2s
#     get_min_results_current = [get_min(λs, θ1, d2_spesific, d3_spesific) for d3_spesific in d3s]
#     min_r_tm_current = [item[1][2] for item in get_min_results_current]
#     min_λ_tm_current = [item[2][2] for item in get_min_results_current]

#     d2_label = round(d2_spesific * 1e9, digits=2)
#     plot!(fig_r_all, d3s .* 1e9, min_r_tm_current, label="d2 = $d2_label nm")
#     plot!(fig_λ_all, d3s .* 1e9, min_λ_tm_current, label="d2 = $d2_label nm")
# end
# display(fig_r_all)
# if save savefig(fig_r_all, fig_path * "min_reflection_all_" * fig_token * ".pdf") end
# display(fig_λ_all)
# if save savefig(fig_λ_all, fig_path * "lambda_min_reflection_all_" * fig_token * ".pdf") end

## Plot plasmon minima and λ for plasmon as a function of d3 for all θds
# fig_r_all = plot(title="Minimum tm reflection for $description at θ = $θd degrees", legend=:topright, xlabel="thickness in nm", ylabel="minimum reflection coefficient", ticks=:native)
# fig_λ_all = plot(title="λ for minimum tm reflection for $description at θ = $θd degrees", legend=:topright, xlabel="thickness in nm", ylabel="λ at minimum reflection coefficient", ticks=:native)
# for θ1_spesific in θs
#     get_min_results_current = [get_min(λs_broadband, θ1_spesific, d2, d3_spesific) for d3_spesific in d3s]
#     min_r_tm_current = [item[1][2] for item in get_min_results_current]
#     min_λ_tm_current = [item[2][2] for item in get_min_results_current]

#     θ_label = round(θ1_spesific * 180 / π, digits=2)
#     plot!(fig_r_all, d3s .* 1e9, min_r_tm_current, label="θ = $θ_label")
#     plot!(fig_λ_all, d3s .* 1e9, min_λ_tm_current, label="θ = $θ_label")
# end
# display(fig_r_all)
# if save savefig(fig_r_all, fig_path * "min_reflection_angles_" * fig_token * ".pdf") end
# display(fig_λ_all)
# if save savefig(fig_λ_all, fig_path * "lambda_min_reflection_angles_" * fig_token * ".pdf") end

## Plot plasmon spectral response for all thicknesses
# fig_all = plot(title="Wavelength respons for all thicknesses for $description at θ = $θd", legend=:outerright, xlabel="wavelength in nm", ylabel="reflection coefficient", ticks=:native)
# for d in d3s
#     systemdata = reduce(vcat, transpose.(system_single.(λs, θ1, d2, d)))
#     d = round(d * 1e9, digits=2)
#     plot!(fig_all, λs .* 1e9, systemdata[:, 1], label="TE, d = $d")
#     plot!(fig_all, λs .* 1e9, systemdata[:, 2], label="TM, d = $d")
# end
# display(fig_all)
# if save savefig(fig_all, fig_path * "responses_" * fig_token * ".pdf") end


println("Finished without errors")