using Plots

using simulator.materials.spesifics
using simulator.targetfigures
using simulator.analyticalmaterials
using simulator.ema


# Simulation parameters start

save = true

n1 = SiO2_core_Sellmeier
n2 = Au_Johnson
n3 = Pd052_unloaded
n4 = Air
n4_ema(λ) = n_Bruggeman(ϵ(n3(λ)), ϵ(Air(λ)), 0.5)
n5 = Air

d1 = 10e-9
d2 = 30e-9
d3 = 10e-9
d4 = 0e-9
d5 = 10e-9

d2smax = 40e-9
d3smax = 40e-9
d4smax = 5e-9
d2ss = 1e-9
d3ss = 1e-10
d4ss = 1e-10
d3smrd2sm = 0.4
d4smrd3sm = 0.1

λ = 1200e-9
λs = [x for x in 600e-9:1e-10:1400e-9]

Hcs = [x for x in 1e-2:1e-2:1e-1]
Pd_c = 0.52

n3s = [Pd034_unloaded, Pd042_unloaded, Pd052_unloaded, Pd073_unloaded, Pd_unloaded]
Pd_cs = [0.34, 0.42, 0.52, 0.73, 1.00]
Pd_d3s = [10e-9, 10e-9, 10e-9, 10e-9, 10e-9]

θd = 45
θds = [x for x in 40:0.001:60]

fig_path = "../rapport/figures/thesisfigures/"

# Simulation parameters end

if save
    gr()
else
    plotly()
end
θ1 = θd * π / 180
θs = θds .* (π / 180)



# # Run thickness simulations
# println("\nScan d3:")
# @time begin
#     system_d3 = make_d3_system(n1, n2, n3, n4, n5, d1, d2, d4, d5, θ1)
#     d3_xs, d3_ysr, d3_ysλ, d3_ysw, d3_ysΔ = scan_plasmon_singleparameter(system_d3, [x for x in 0:d3ss:min(d3smax, d2*d3smrd2sm)], λs, 2)
# end
# d3_xsnm = d3_xs *= 1e9
# fig_d3r = plot(d3_xsnm, d3_ysr, title="Plasmon reflection minimum", xaxis="Palladium thickness [nm]", yaxis="Power reflection coefficient minimum", legend=false, tiks=:native)
# fig_d3λ = plot(d3_xsnm, d3_ysλ .* 1e9, title="Plasmon minimum wavelength", xaxis="Palladium thickness [nm]", yaxis="Dip wavelength [nm]", legend=false, tiks=:native)
# fig_d3w = plot(d3_xsnm, d3_ysw .* 1e9, title="Plasmon dip full width half maximum", xaxis="Palladium thickness [nm]", yaxis="Dip full width half maximum [nm]", legend=false, tiks=:native)
# fig_d3Δ = plot(d3_xsnm, d3_ysΔ, title="Plasmon dip peak to peak", xaxis="Palladium thickness [nm]", yaxis="Dip peak to peak", legend=false, tiks=:native)
# if save
#     savefig(fig_d3r, fig_path * "d3r.pdf")
#     savefig(fig_d3λ, fig_path * "d3l.pdf")
#     savefig(fig_d3w, fig_path * "d3w.pdf")
#     savefig(fig_d3Δ, fig_path * "d3d.pdf")
# else
#     display(fig_d3r)
#     display(fig_d3λ)
#     display(fig_d3w)
#     display(fig_d3Δ)
# end


# println("\nScan d2 and d3:")
# @time begin
#     system_d23 = make_d2_dualparameter_system(n1, n2, n3, n4, n5, d1, d4, d5, θ1)
#     d23_xs, d23_ysr, d23_ysλ, d23_ysw, d23_ysΔ = scan_plasmon_dualthicknesses(system_d23, 0, d2smax, d2ss, 0, d3smax, d2ss, λs, 1, d3smrd2sm)
# end
# d23_xsnm = d23_xs .* 1e9
# fig_d23r = plot(d23_xsnm, d23_ysr .* 1e9, title="Minimum plasmon dip", xaxis="Gold thickness [nm]", yaxis="Palladium thickness [nm]", legend=false, tiks=:native)
# fig_d23w = plot(d23_xsnm, d23_ysw .* 1e9, title="Minimum dip width", xaxis="Gold thickness [nm]", yaxis="Palladium thickness [nm]", legend=false, tiks=:native)
# fig_d23Δ = plot(d23_xsnm, d23_ysΔ .* 1e9, title="Maximum dip peak to peak", xaxis="Gold thickness [nm]", yaxis="Palladium thickness [nm]", legend=false, tiks=:native)
# if save
#     savefig(fig_d23r, fig_path * "d23r.pdf")
#     savefig(fig_d23w, fig_path * "d23w.pdf")
#     savefig(fig_d23Δ, fig_path * "d23d.pdf")
# else
#     display(fig_d23r)
#     display(fig_d23w)
#     display(fig_d23Δ)
# end



# # Run demo of structure with and without gold
# println("\nSingle plasmon with without gold")
# @time begin
#     system_with = make_layered_tm_system([n1, n2, n3, n4, n5], [d1, 30e-9, 10e-9, d4, d5], θ1)
#     system_without = make_layered_tm_system([n1, n2, n3, n4, n5], [d1, 0e-9, 10e-9, d4, d5], θ1)
#     with_xs, with_ys = scann_singleparameter(system_with, λs)
#     without_xs, without_ys = scann_singleparameter(system_without, λs)

#     system_with_θ = make_layered_θ_tm_system([n1, n2, n3, n4, n5], [d1, 30e-9, 10e-9, d4, d5], λ)
#     system_without_θ = make_layered_θ_tm_system([n1, n2, n3, n4, n5], [d1, 0e-9, 10e-9, d4, d5], λ)
#     with_xsθ, with_ysθ = scann_singleparameter(system_with_θ, θs)
#     without_xsθ, without_ysθ = scann_singleparameter(system_without_θ, θs)
# end
# fig_with = plot(with_xs .* 1e9, with_ys, title="Plasmon structure with and without gold", xaxis="Incident wavelength [nm]", yaxis="Reflection coefficient", label="10 nm palladium on 30 nm gold", legend=:bottomright, tiks=:native)
# plot!(fig_with, without_xs .* 1e9, without_ys, label="10 nm palladium")
# fig_with_θ = plot(θds, with_ysθ, title="Plasmon structure with and without gold", xaxis="Incident angle [degrees]", yaxis="Reflection coefficient", label="10 nm palladium on 30 nm gold", legend=:bottomright, tiks=:native)
# plot!(fig_with_θ, θds, without_ysθ, label="10 nm palladium")
# if save
#     savefig(fig_with, fig_path * "withwithout.pdf")
#     savefig(fig_with_θ, fig_path * "withwithouttheta.pdf")
# else
#     display(fig_with)
#     display(fig_with_θ)
# end



# # Run hydrogen simulation
# println("\nScan hydrogen:")
# @time begin
#     system_Hc = make_hc_system(n1, n2, n3, n4, n5, d1, d2, d3, d4, d5, θ1, Pd_c)
#     Hc_xs, Hc_ysr, Hc_ysλ, Hc_ysw, Hc_ysΔ = scan_plasmon_singleparameter(system_Hc, Hcs, λs, 2)
# end
# Hc_xsp = Hc_xs .* 1e2
# fig_Hcr = plot(Hc_xsp, Hc_ysr, title="Plasmon reflection minimum", xaxis="Hydrogen concentration [%]", yaxis="Power reflection coefficient minimum", legend=false, tiks=:native)
# fig_Hcλ = plot(Hc_xsp, Hc_ysλ .* 1e9, title="Plasmon minimum wavelength", xaxis="Hydrogen concentration [%]", yaxis="Dip wavelength [nm]", legend=false, tiks=:native)
# fig_Hcw = plot(Hc_xsp, Hc_ysw .* 1e9, title="Plasmon dip full width half maximum", xaxis="Hydrogen concentration [%]", yaxis="Dip full width half maximum [nm]", legend=false, tiks=:native)
# fig_HcΔ = plot(Hc_xsp, Hc_ysΔ, title="Plasmon dip peak to peak", xaxis="Hydrogen concentration [%]", yaxis="Dip peak to peak", legend=false, tiks=:native)
# if save
#     savefig(fig_Hcr, fig_path * "Hcr.pdf")
#     savefig(fig_Hcλ, fig_path * "Hcl.pdf")
#     savefig(fig_Hcw, fig_path * "Hcw.pdf")
#     savefig(fig_HcΔ, fig_path * "Hcd.pdf")
# else
#     display(fig_Hcr)
#     display(fig_Hcλ)
#     display(fig_Hcw)
#     display(fig_HcΔ)
# end



# # Run alloy comparisons
# println("\nScan alloys d3:")
# fig_alloys_d3r = plot(title="Plasmon reflection minimum", xaxis="Palladium thickness [nm]", yaxis="Power reflection coefficient minimum", legend=:bottomright, tiks=:native)
# fig_alloys_d3λ = plot(title="Plasmon minimum wavelength", xaxis="Palladium thickness [nm]", yaxis="Dip wavelength [nm]", legend=:bottomleft, tiks=:native)
# fig_alloys_d3w = plot(title="Plasmon dip full width half maximum", xaxis="Palladium thickness [nm]", yaxis="Dip full width half maximum [nm]", legend=:bottomright, tiks=:native)
# fig_alloys_d3Δ = plot(title="Plasmon dip peak to peak", xaxis="Palladium thickness [nm]", yaxis="Dip peak to peak", legend=:topright, tiks=:native)
# @time begin
#     for (n3_local, Pd_c_local) in zip(n3s, Pd_cs)
#         @show Pd_c_local
#         system_alloysd3 = make_d3_system(n1, n2, n3_local, n4, n5, d1, d2, d4, d5, θ1)
#         alloysd3_xs, alloysd3_ysr, alloysd3_ysλ, alloysd3_ysw, alloysd3_ysΔ = scan_plasmon_singleparameter(system_alloysd3, [x for x in 0:d3ss:min(d3smax, d2*d3smrd2sm)], λs, 2)
        
#         alloysd3_xsnm = alloysd3_xs .* 1e9
#         plot!(fig_alloys_d3r, alloysd3_xsnm, alloysd3_ysr, label=(string(Pd_c_local) * "% palladium"))
#         plot!(fig_alloys_d3λ, alloysd3_xsnm, alloysd3_ysλ .* 1e9 , label=(string(Pd_c_local) * "% palladium"))
#         plot!(fig_alloys_d3w, alloysd3_xsnm, alloysd3_ysw .* 1e9 , label=(string(Pd_c_local) * "% palladium"))
#         plot!(fig_alloys_d3Δ, alloysd3_xsnm, alloysd3_ysΔ, label=(string(Pd_c_local) * "% palladium"))
#     end
# end
# if save
#     savefig(fig_alloys_d3r, fig_path * "alloys_d3r.pdf")
#     savefig(fig_alloys_d3λ, fig_path * "alloys_d3l.pdf")
#     savefig(fig_alloys_d3w, fig_path * "alloys_d3w.pdf")
#     savefig(fig_alloys_d3Δ, fig_path * "alloys_d3d.pdf")
# else
#     display(fig_alloys_d3r)
#     display(fig_alloys_d3λ)
#     display(fig_alloys_d3w)
#     display(fig_alloys_d3Δ)
# end


# println("\nScan alloys hydrogen:")
# fig_alloys_Hcr = plot(title="Plasmon reflection minimum", xaxis="Hydrogen concentration [%]", yaxis="Power reflection coefficient minimum", legend=:bottomright, tiks=:native)
# fig_alloys_Hcλ = plot(title="Plasmon minimum wavelength", xaxis="Hydrogen concentration [%]", yaxis="Dip wavelength [nm]", legend=:bottomright, tiks=:native)
# fig_alloys_Hcw = plot(title="Plasmon dip full width half maximum", xaxis="Hydrogen concentration [%]", yaxis="Dip full width half maximum [nm]", legend=:bottomright, tiks=:native)
# fig_alloys_HcΔ = plot(title="Plasmon dip peak to peak", xaxis="Hydrogen concentration [%]", yaxis="Dip peak to peak", legend=:bottomright, tiks=:native)
# @time begin
#     for (n3_local, Pd_c_local, d3_local) in collect(zip(n3s, Pd_cs, Pd_d3s))[2:length(n3s)]
#         @show Pd_c_local
#         system_alloys_Hc = make_hc_system(n1, n2, n3_local, n4, n5, d1, d2, d3_local, d4, d5, θ1, Pd_c_local)
#         alloysHc_xs, alloysHc_ysr, alloysHc_ysλ, alloysHc_ysw, alloysHc_ysΔ = scan_plasmon_singleparameter(system_alloys_Hc, Hcs, λs, 2)
        
#         alloysHc_xsp = alloysHc_xs .* 1e2
#         plot!(fig_alloys_Hcr, alloysHc_xsp, alloysHc_ysr, label=(string(Pd_c_local) * "% palladium"))
#         plot!(fig_alloys_Hcλ, alloysHc_xsp, alloysHc_ysλ .* 1e9 , label=(string(Pd_c_local) * "% palladium"))
#         plot!(fig_alloys_Hcw, alloysHc_xsp, alloysHc_ysw .* 1e9 , label=(string(Pd_c_local) * "% palladium"))
#         plot!(fig_alloys_HcΔ, alloysHc_xsp, alloysHc_ysΔ, label=(string(Pd_c_local) * "% palladium"))
#     end
# end
# if save
#     savefig(fig_alloys_Hcr, fig_path * "alloys_Hcr.pdf")
#     savefig(fig_alloys_Hcλ, fig_path * "alloys_Hcl.pdf")
#     savefig(fig_alloys_Hcw, fig_path * "alloys_Hcw.pdf")
#     savefig(fig_alloys_HcΔ, fig_path * "alloys_Hcd.pdf")
# else
#     display(fig_alloys_Hcr)
#     display(fig_alloys_Hcλ)
#     display(fig_alloys_Hcw)
#     display(fig_alloys_HcΔ)
# end



# # Run roughness simulations
# println("\nScan d4 roughness:")
# @time begin
#     system_d4_ema = make_d4_system(n1, n2, n3, n4_ema, n5, d1, d2, d3, d5, θ1)
#     d4_ema_xs, d4_ema_ysr, d4_ema_ysλ, d4_ema_ysw, d4_ema_ysΔ = scan_plasmon_singleparameter(system_d4_ema, [x for x in 0:d4ss:min(d4smax, d3*d4smrd3sm)], λs, 2)
# end
# d4_ema_xsnm = d4_ema_xs *= 1e9
# fig_d4_emar = plot(d4_ema_xsnm, d4_ema_ysr, title="Plasmon reflection minimum", xaxis="Roughness ema thickness [nm]", yaxis="Power reflection coefficient minimum", legend=false, tiks=:native)
# fig_d4_emaλ = plot(d4_ema_xsnm, d4_ema_ysλ .* 1e9, title="Plasmon minimum wavelength", xaxis="Roughness ema thickness [nm]", yaxis="Dip wavelength [nm]", legend=false, tiks=:native)
# fig_d4_emaw = plot(d4_ema_xsnm, d4_ema_ysw .* 1e9, title="Plasmon dip full width half maximum", xaxis="Roughness ema thickness [nm]", yaxis="Dip full width half maximum [nm]", legend=false, tiks=:native)
# fig_d4_emaΔ = plot(d4_ema_xsnm, d4_ema_ysΔ, title="Plasmon dip peak to peak", xaxis="Roughness ema thickness [nm]", yaxis="Dip peak to peak", legend=false, tiks=:native)
# if save
#     savefig(fig_d4_emar, fig_path * "d4_emar.pdf")
#     savefig(fig_d4_emaλ, fig_path * "d4_emal.pdf")
#     savefig(fig_d4_emaw, fig_path * "d4_emaw.pdf")
#     savefig(fig_d4_emaΔ, fig_path * "d4_emad.pdf")
# else
#     display(fig_d4_emar)
#     display(fig_d4_emaλ)
#     display(fig_d4_emaw)
#     display(fig_d4_emaΔ)
# end
