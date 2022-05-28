using Plots
using simulator.materials.spesifics

# plotly()
gr()

λs = 500e-9:1e-9:1600e-9


fig_n_Au = plot(title="Real part of gold refractive indicies", legend=:outerright, xlabel="Wavelength [nm]", ylabel="n", ticks=:native)
plot!(fig_n_Au, λs .* 1e9, real(Au_Johnson.(λs)), label="Johnson")
plot!(fig_n_Au, λs .* 1e9, real(Au_Werner.(λs)), label="Werner")
plot!(fig_n_Au, λs .* 1e9, real(Au_unloaded.(λs)), label="Palm 2019")
plot!(fig_n_Au, λs .* 1e9, real(Au_McPeak.(λs)), label="McPeak")
plot!(fig_n_Au, λs .* 1e9, real(Au_Babar.(λs)), label="Babar")
plot!(fig_n_Au, λs .* 1e9, real(Au_OlmonEvaporated.(λs)), label="Olmon evaporated")
plot!(fig_n_Au, λs .* 1e9, real(Au_OlmonSingleChrystaline.(λs)), label="Olmon single chrystaline")
plot!(fig_n_Au, λs .* 1e9, real(Au_OlmonTemplateStripped.(λs)), label="Olmon template stripped")
plot!(fig_n_Au, λs .* 1e9, real(Au_11nm_Rosenblatt.(λs)), label="Rosenblatt 11 nm film")
plot!(fig_n_Au, λs .* 1e9, real(Au_21nm_Rosenblatt.(λs)), label="Rosenblatt 21 nm film")
plot!(fig_n_Au, λs .* 1e9, real(Au_44nm_Rosenblatt.(λs)), label="Rosenblatt 44 nm film")
plot!(fig_n_Au, λs .* 1e9, real(SiO2_core_Sellmeier.(λs)), label="SiO2 core")
display(fig_n_Au)
# savefig(fig_n_Au, "../rapport/figures/IndexSourceComparisons/Au_n.pdf")

fig_k_Au = plot(title="Imaginary part of gold refractive indicies", legend=:bottomleft, xlabel="Wavelength [nm]", ylabel="k", ticks=:native)
plot!(fig_k_Au, λs .* 1e9, imag(Au_Johnson.(λs)), label="Johnson")
plot!(fig_k_Au, λs .* 1e9, imag(Au_Werner.(λs)), label="Werner")
plot!(fig_k_Au, λs .* 1e9, imag(Au_unloaded.(λs)), label="Palm 2019")
plot!(fig_k_Au, λs .* 1e9, imag(Au_McPeak.(λs)), label="McPeak")
plot!(fig_k_Au, λs .* 1e9, imag(Au_Babar.(λs)), label="Babar")
plot!(fig_k_Au, λs .* 1e9, imag(Au_OlmonEvaporated.(λs)), label="Olmon evaporated")
plot!(fig_k_Au, λs .* 1e9, imag(Au_OlmonSingleChrystaline.(λs)), label="Olmon single chrystaline")
plot!(fig_k_Au, λs .* 1e9, imag(Au_OlmonTemplateStripped.(λs)), label="Olmon template stripped")
plot!(fig_k_Au, λs .* 1e9, imag(Au_11nm_Rosenblatt.(λs)), label="Rosenblatt 11 nm film")
plot!(fig_k_Au, λs .* 1e9, imag(Au_21nm_Rosenblatt.(λs)), label="Rosenblatt 21 nm film")
plot!(fig_k_Au, λs .* 1e9, imag(Au_44nm_Rosenblatt.(λs)), label="Rosenblatt 44 nm film")
plot!(fig_k_Au, λs .* 1e9, imag(SiO2_core_Sellmeier.(λs)), label="SiO2 core")
display(fig_k_Au)
# savefig(fig_k_Au, "../rapport/figures/IndexSourceComparisons/Au_k.pdf")


fig_n_Pd = plot(title="Real part of palladium refractive indicies", legend=:outerright, xlabel="Wavelength [nm]", ylabel="n", ticks=:native)
plot!(fig_n_Pd, λs .* 1e9, real(Pd_Johnson.(λs)), label="Johnson")
plot!(fig_n_Pd, λs .* 1e9, real(Pd_Werner.(λs)), label="Werner")
plot!(fig_n_Pd, λs .* 1e9, real(Pd_Palm_2018.(λs)), label="Palm 2018")
plot!(fig_n_Pd, λs .* 1e9, real(Pd_unloaded.(λs)), label="Palm 2019 unloaded")
plot!(fig_n_Pd, λs .* 1e9, real(Pd_loaded.(λs)), label="Palm 2019 loaded")
plot!(fig_n_Pd, λs .* 1e9, real(SiO2_core_Sellmeier.(λs)), label="SiO2 core")
display(fig_n_Pd)
# savefig(fig_n_Pd, "../rapport/figures/IndexSourceComparisons/Pd_n.pdf")

fig_k_Pd = plot(title="Imaginary part of palladium refractive indicies", legend=:bottomleft, xlabel="Wavelength [nm]", ylabel="k", ticks=:native)
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_Johnson.(λs)), label="Johnson")
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_Werner.(λs)), label="Werner")
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_Palm_2018.(λs)), label="Palm 2018")
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_unloaded.(λs)), label="Palm 2019 unloaded")
plot!(fig_k_Pd, λs .* 1e9, imag(Pd_loaded.(λs)), label="Palm 2019 loaded")
plot!(fig_k_Pd, λs .* 1e9, imag(SiO2_core_Sellmeier.(λs)), label="SiO2 core")
display(fig_k_Pd)
# savefig(fig_k_Pd, "../rapport/figures/IndexSourceComparisons/Pd_k.pdf")


fig_n_lab = plot(title="Real part of refractive index for samples", legend=:outerright, xlabel="Wavelength [nm]", ylabel="n", ticks=:native)
plot!(fig_n_lab, λs .* 1e9, real(Au_Johnson.(λs)), label="Au Johnson")
plot!(fig_n_lab, λs .* 1e9, real(Pd_Johnson.(λs)), label="Pd Johnson")

plot!(fig_n_lab, λs .* 1e9, real(Pd000S5.(λs)), label="S5 0% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd100S20.(λs)), label="S20 100% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd070S9to17.(λs)), label="S9-17 70% Pd")

plot!(fig_n_lab, λs .* 1e9, real(Pd041S21.(λs)), label="S21 41% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd051S22.(λs)), label="S22 51% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd061S23.(λs)), label="S23 61% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd069S24.(λs)), label="S24 69% Pd")

plot!(fig_n_lab, λs .* 1e9, real(Pd064S7and18.(λs)), label="S7, 18 64% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd053S6.(λs)), label="S6 53% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd074S8drude.(λs)), label="S8 drude 74% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd061S19.(λs)), label="S19 61% Pd")
plot!(fig_n_lab, λs .* 1e9, real(Pd061S19drude.(λs)), label="S19 drude 61% Pd")
display(fig_n_lab)
# savefig(fig_n_lab, "../rapport/figures/thesisfigures/Samples_n.pdf")

fig_k_lab = plot(title="Imaginary part of refractive index for samples", legend=:bottomleft, xlabel="Wavelength [nm]", ylabel="k", ticks=:native)
plot!(fig_k_lab, λs .* 1e9, imag(Au_Johnson.(λs)), label="Au Johnson")
plot!(fig_k_lab, λs .* 1e9, imag(Pd_Johnson.(λs)), label="Pd Johnson")

plot!(fig_k_lab, λs .* 1e9, imag(Pd000S5.(λs)), label="S5 0% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd100S20.(λs)), label="S20 100% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd070S9to17.(λs)), label="S9-17 70% Pd")

plot!(fig_k_lab, λs .* 1e9, imag(Pd041S21.(λs)), label="S21 41% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd051S22.(λs)), label="S22 51% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd061S23.(λs)), label="S23 61% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd069S24.(λs)), label="S24 69% Pd")

plot!(fig_k_lab, λs .* 1e9, imag(Pd064S7and18.(λs)), label="S7, 18 64% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd053S6.(λs)), label="S6 53% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd074S8drude.(λs)), label="S8 drude 74% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd061S19.(λs)), label="S19 61% Pd")
plot!(fig_k_lab, λs .* 1e9, imag(Pd061S19drude.(λs)), label="S19 drude 61% Pd")
display(fig_k_lab)
# savefig(fig_k_lab, "../rapport/figures/thesisfigures/Samples_k.pdf")


fig_n_comp = plot(title="Real part of refractive index for samples", legend=:outerright, xlabel="Wavelength [nm]", ylabel="n", ticks=:native)
plot!(fig_n_comp, λs .* 1e9, real(SiO2_core_Sellmeier.(λs)), label="SiO2 core")

plot!(fig_n_comp, λs .* 1e9, real(Pd041S21.(λs)), label="S21 41% Pd")
plot!(fig_n_comp, λs .* 1e9, real(Pd051S22.(λs)), label="S22 51% Pd")
plot!(fig_n_comp, λs .* 1e9, real(Pd053S6.(λs)), label="S6 53% Pd")
plot!(fig_n_comp, λs .* 1e9, real(Pd070S9to17.(λs)), label="S9-17 70% Pd")
plot!(fig_n_comp, λs .* 1e9, real(Pd074S8drude.(λs)), label="S8 drude 74% Pd")

plot!(fig_n_comp, λs .* 1e9, real(Pd042_unloaded.(λs)), label="Palm 42% Pd")
plot!(fig_n_comp, λs .* 1e9, real(Pd052_unloaded.(λs)), label="Palm 52% Pd")
plot!(fig_n_comp, λs .* 1e9, real(Pd073_unloaded.(λs)), label="Palm 73% Pd")
display(fig_n_comp)
# savefig(fig_n_comp, "../rapport/figures/thesisfigures/Samples_n_comp.pdf")

fig_k_comp = plot(title="Imaginary part of refractive index for samples", legend=:bottomleft, xlabel="Wavelength [nm]", ylabel="k", ticks=:native)
plot!(fig_k_comp, λs .* 1e9, imag(SiO2_core_Sellmeier.(λs)), label="SiO2 core")

plot!(fig_k_comp, λs .* 1e9, imag(Pd041S21.(λs)), label="S21 41% Pd")
plot!(fig_k_comp, λs .* 1e9, imag(Pd051S22.(λs)), label="S22 51% Pd")
plot!(fig_k_comp, λs .* 1e9, imag(Pd053S6.(λs)), label="S6 53% Pd")
plot!(fig_k_comp, λs .* 1e9, imag(Pd070S9to17.(λs)), label="S9-17 70% Pd")
plot!(fig_k_comp, λs .* 1e9, imag(Pd074S8drude.(λs)), label="S8 drude 74% Pd")

plot!(fig_k_comp, λs .* 1e9, imag(Pd042_unloaded.(λs)), label="Palm 42% Pd")
plot!(fig_k_comp, λs .* 1e9, imag(Pd052_unloaded.(λs)), label="Palm 52% Pd")
plot!(fig_k_comp, λs .* 1e9, imag(Pd073_unloaded.(λs)), label="Palm 73% Pd")
display(fig_k_comp)
# savefig(fig_k_comp, "../rapport/figures/thesisfigures/Samples_k_comp.pdf")
