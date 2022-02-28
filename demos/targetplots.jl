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
d2s = [x for x in 0:1e-9:40e-9]
d3s = [x for x in 0:1e-10:40e-9]

λ = 1200e-9
λs = [x for x in 600e-9:1e-9:1400e-9]

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
xs, ys = scann_singleparameter(system_single, λs)
end
fig_single = plot(xs .* 1e9, ys, title="$description, Single scan spectra", xaxis="λ [nm]", yaxis="Power reflection coefficient", legend=false)
display(fig_single)
if save savefig(fig_single, fig_path * "single_" * fig_token * ".pdf") end

# Plot d3 scan:
println("\nScan d3 timing:")
@time begin
system_d3 = make_d3_system(n1, n2, n3, n4, n5, d1, d2, d4, d5, θ1)
xs, ysr, _ = scan_plasmon_singleparameter(system_d3, d3s, λs)
end
fig_d3 = plot(xs .* 1e9, ysr, title="$description, Active layer scan", xaxis="d [nm]", yaxis="Power reflection coefficient minima", legend=false)
display(fig_d3)
if save savefig(fig_d3, fig_path * "d3_" * fig_token * ".pdf") end

# Plot d2 on d3 scan:
println("\nScan d2 and d3 timing:")
@time begin
system_d2 = make_d2_dualparameter_system(n1, n2, n3, n4, n5, d1, d4, d5, θ1)
xs_d2, ys_d3, _ = scan_plasmon_dualparameter(system_d2, d2s, d3s, λs)
end
fig_d2 = plot(xs_d2 .* 1e9, ys_d3 .* 1e9, title="$description, Double layer scan", xaxis="Au thickness [nm]", yaxis="Pd thickness [nm]", legend=false)
display(fig_d2)
if save savefig(fig_d2, fig_path * "d2-d3_" * fig_token * ".pdf") end
