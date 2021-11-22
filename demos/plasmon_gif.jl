using Plots

using simulator.materials.spesifics
using simulator.fresnelltools
using simulator.matrixcore


λ = 1200e-9

n1 = 1.5
nm = Pd(λ)
n2 = 1.2

θs = 0.7:1e-4:1


function point(θ, n2, which)
    i_1_2_te, i_1_2_tm, θm = FresnellBoundrary(n1, nm, θ)
    bulk = FresnellSlab(nm, 2*π/λ, 50e-9, θm)
    i_2_3_te, i_2_3_tm, _ = FresnellBoundrary(nm, n2, θm)

    system_te = CascadeScattering([i_1_2_te, bulk, i_2_3_te])
    system_tm = CascadeScattering([i_1_2_tm, bulk, i_2_3_tm])

    te_r = abs(system_te[2, 1])^2
    tm_r = abs(system_tm[2, 1])^2

    [te_r, tm_r][which]
end

anim = @animate for n ∈ 1:0.001:1.2
    plot(θs, point.(θs, n, 1), ylim = [0.5, 1], title="Plasmon reflectance in varying index environment", legend=:bottomright, label="TE")
    plot!(θs, point.(θs, n, 2), label="TM")
end

mp4(anim, "~/Downloads/temp.mp4", fps=30)