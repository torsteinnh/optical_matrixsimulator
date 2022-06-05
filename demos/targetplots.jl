using Plots

using simulator.materials.spesifics
using simulator.targetfigures


# Simulation parameters start

save = false

n1 = SiO2_core_Sellmeier
n2 = Au_Johnson
n3 = Pd052_unloaded
n4 = Air
n5 = Air


d1 = 10e-9
d2 = 30e-9
d3 = 5e-9
d4 = 10e-9
d5 = 10e-9
d2smax = 40e-9
d3smax = 40e-9
d2ss = 1e-9
d3ss = 1e-10
d3smrd2sm = 0.4

λ = 1200e-9
λs = [x for x in 600e-9:1e-10:1400e-9]

Hcs = [x for x in 1e-2:1e-2:1e-1]
Pd_c = 0.52

θd = 45

description = "52 unloaded Pd on Au"
fig_path = "../rapport/figures/targetplots/"

# Simulation parameters end


fig_token = replace(description, " " => "_")

if save
    gr()
else
    plotly()
end
θ1 = θd * π / 180


# Plot singe line:
println("\nSingle line timing:")
@time begin
system_single = make_layered_tm_system([n1, n2, n3, n4, n5], [d1, d2, d3, d4, d5], θ1)
xs, ys = scan_singleparameter(system_single, λs)
end
fig_single = plot(xs .* 1e9, ys, title="$description, Single scan spectra", xaxis="λ [nm]", yaxis="Power reflection coefficient", legend=false)
display(fig_single)
if save savefig(fig_single, fig_path * "single_" * fig_token * ".pdf") end

# Plot hydrogen concentration scan:
println("\nScan hydrogen timing:")
@time begin
system_Hc = make_hc_system(n1, n2, n3, n4, n5, d1, d2, d3, d4, d5, θ1, Pd_c)
xs, ysr, ysλ, ysw, ysΔ = scan_plasmon_singleparameter(system_Hc, Hcs, λs, 2)
end
fig_d3r = plot(xs .* 1e2, ysr, title="$description, hydrogen scan", xaxis="Hydrogen concentration [%]", yaxis="Power reflection coefficient minima", legend=false)
fig_d3λ = plot(xs .* 1e2, ysλ .* 1e9, title="$description, hydrogen scan", xaxis="Hydrogen concentration [%]", yaxis="λ [nm]", legend=false)
fig_d3w = plot(xs .* 1e2, ysw .* 1e9, title="$description, hydrogen scan", xaxis="Hydrogen concentration [%]", yaxis="λ plasmon dip FWHM [nm]", legend=false)
fig_d3Δ = plot(xs .* 1e2, ysΔ, title="$description, hydrogen scan", xaxis="Hydrogen concentration [%]", yaxis="Plasmon peak power reflection Δ", legend=false)
display(fig_d3r)
if save savefig(fig_d3r, fig_path * "Hcr_" * fig_token * ".pdf") end
display(fig_d3λ)
if save savefig(fig_d3λ, fig_path * "Hcl_" * fig_token * ".pdf") end
display(fig_d3w)
if save savefig(fig_d3w, fig_path * "Hcw_" * fig_token * ".pdf") end
display(fig_d3Δ)
if save savefig(fig_d3Δ, fig_path * "Hcd_" * fig_token * ".pdf") end

# Plot d3 scan:
println("\nScan d3 timing:")
@time begin
system_d3 = make_d3_system(n1, n2, n3, n4, n5, d1, d2, d4, d5, θ1)
xs, ysr, ysλ, ysw, ysΔ = scan_plasmon_singleparameter(system_d3, [x for x in 0:d3ss:(d2*d3smrd2sm)], λs, 2)
end
fig_d3r = plot(xs .* 1e9, ysr, title="$description, Active layer scan", xaxis="d [nm]", yaxis="Power reflection coefficient minima", legend=false)
fig_d3λ = plot(xs .* 1e9, ysλ .* 1e9, title="$description, Active layer scan", xaxis="d [nm]", yaxis="λ [nm]", legend=false)
fig_d3w = plot(xs .* 1e9, ysw .* 1e9, title="$description, Active layer scan", xaxis="d [nm]", yaxis="λ plasmon dip FWHM [nm]", legend=false)
fig_d3Δ = plot(xs .* 1e9, ysΔ, title="$description, Active layer scan", xaxis="d [nm]", yaxis="Plasmon peak power reflection Δ", legend=false)
display(fig_d3r)
if save savefig(fig_d3r, fig_path * "d3r_" * fig_token * ".pdf") end
display(fig_d3λ)
if save savefig(fig_d3λ, fig_path * "d3l_" * fig_token * ".pdf") end
display(fig_d3w)
if save savefig(fig_d3w, fig_path * "d3w_" * fig_token * ".pdf") end
display(fig_d3Δ)
if save savefig(fig_d3Δ, fig_path * "d3d_" * fig_token * ".pdf") end

# Plot d2 on d3 scan:
println("\nScan d2 and d3 timing:")
@time begin
system_d2 = make_d2_dualparameter_system(n1, n2, n3, n4, n5, d1, d4, d5, θ1)
xs_d2, ys_d3R, _, ys_d3ν, ys_d3Δ = scan_plasmon_dualthicknesses(system_d2, 0, d2smax, d2ss, 0, d3smax, d3ss, λs, 1, d3smrd2sm)
end
fig_d2R = plot(xs_d2 .* 1e9, ys_d3R .* 1e9, title="$description, minimum plasmon dip", xaxis="Au thickness [nm]", yaxis="Pd thickness [nm]", legend=false)
display(fig_d2R)
if save savefig(fig_d2R, fig_path * "d2-d3-R_" * fig_token * ".pdf") end
fig_d2ν = plot(xs_d2 .* 1e9, ys_d3ν .* 1e9, title="$description, minimum plasmon FWHM", xaxis="Au thickness [nm]", yaxis="Pd thickness [nm]", legend=false)
display(fig_d2ν)
if save savefig(fig_d2ν, fig_path * "d2-d3-FWHM_" * fig_token * ".pdf") end
fig_d2R = plot(xs_d2 .* 1e9, ys_d3Δ .* 1e9, title="$description, maximum plasmon peak to peak", xaxis="Au thickness [nm]", yaxis="Pd thickness [nm]", legend=false)
display(fig_d2R)
if save savefig(fig_d2R, fig_path * "d2-d3-P2P_" * fig_token * ".pdf") end
