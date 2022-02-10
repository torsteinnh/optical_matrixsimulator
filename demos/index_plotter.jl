using Plots
using simulator.materials.spesifics

gr()

λs = 400e-9:1e-9:1500e-9


fig_n_Au = plot(title="Real part of gold refractive indicies", legend=:outerright, xlabel="Wavelength [nm]", ylabel="n", ticks=:native)
plot!(fig_n_Au, λs .* 1e9, real(Au_Johnson.(λs)), label="Johnson")
plot!(fig_n_Au, λs .* 1e9, real(Au_Werner.(λs)), label="Werner")
plot!(fig_n_Au, λs .* 1e9, real(Au_unloaded.(λs)), label="Palm 2019")
plot!(fig_n_Au, λs .* 1e9, real(SiO2_core_Sellmeier.(λs)), label="SiO2 core")
display(fig_n_Au)
savefig(fig_n_Au, "../rapport/figures/IndexSourceComparisons/Au_n.pdf")

fig_k_Au = plot(title="Real part of gold refractive indicies", legend=:outerright, xlabel="Wavelength [nm]", ylabel="k", ticks=:native)
plot!(fig_k_Au, λs .* 1e9, imag(Au_Johnson.(λs)), label="Johnson")
plot!(fig_k_Au, λs .* 1e9, imag(Au_Werner.(λs)), label="Werner")
plot!(fig_k_Au, λs .* 1e9, imag(Au_unloaded.(λs)), label="Palm 2019")
plot!(fig_k_Au, λs .* 1e9, imag(SiO2_core_Sellmeier.(λs)), label="SiO2 core")
display(fig_k_Au)
savefig(fig_k_Au, "../rapport/figures/IndexSourceComparisons/Au_k.pdf")


fig_n_Pd = plot(title="Real part of palladium refractive indicies", legend=:outerright, xlabel="Wavelength [nm]", ylabel="n", ticks=:native)
plot!(fig_n_Pd, λs .* 1e9, real(Pd_Johnson.(λs)), label="Johnson")
plot!(fig_n_Pd, λs .* 1e9, real(Pd_Werner.(λs)), label="Werner")
plot!(fig_n_Pd, λs .* 1e9, real(Pd_Palm_2018.(λs)), label="Palm 2018")
plot!(fig_n_Pd, λs .* 1e9, real(Pd_unloaded.(λs)), label="Palm 2019 unloaded")
plot!(fig_n_Pd, λs .* 1e9, real(Pd_loaded.(λs)), label="Palm 2019 loaded")
plot!(fig_n_Pd, λs .* 1e9, real(SiO2_core_Sellmeier.(λs)), label="SiO2 core")
display(fig_n_Pd)
savefig(fig_n_Pd, "../rapport/figures/IndexSourceComparisons/Pd_n.pdf")

fig_k_Pd = plot(title="Real part of palladium refractive indicies", legend=:outerright, xlabel="Wavelength [nm]", ylabel="k", ticks=:native)
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_Johnson.(λs)), label="Johnson")
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_Werner.(λs)), label="Werner")
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_Palm_2018.(λs)), label="Palm 2018")
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_unloaded.(λs)), label="Palm 2019 unloaded")
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_loaded.(λs)), label="Palm 2019 loaded")
plot!(fig_k_Pd, λs .* 1e9, imag(SiO2_core_Sellmeier.(λs)), label="SiO2 core")
display(fig_k_Pd)
savefig(fig_k_Pd, "../rapport/figures/IndexSourceComparisons/Pd_k.pdf")
