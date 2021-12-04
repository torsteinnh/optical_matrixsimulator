using Plots

using simulator.fresnelltools
using simulator.matrixcore
using simulator.materials.spesifics

save = false

n_1 = l -> 1.5
n_2 = Au
n_3 = l -> 1
n_4 = l -> 1
n_5 = l -> 1

d_1 = 30e-9
d_2 = 50e-9
d_3 = 10e-9
d_4 = 10e-9
d_5 = 10e-9

λ = 1200e-9
θd = 42.23

step = 1e-10
expression_total(Up, Un) = abs(Up + Un)^2
expression_forward(Up, _) = abs(Up)^2
expression_backward(_, Un) = abs(Un)^2
expression_borders(u) = abs(u)^2

description = "silver silica 4"
θds = 0:1e-2:90
λs = 400e-9:1e-9:1800e-9

fig_path = "../rapport/figures/interrogator_example/"
fig_token = replace(description, " " => "_")



if save
    gr()
else
    plotly()
end
θ1 = θd * π/180
θs = θds .* (π/180)

function system_slabs(θ1, λ, select)
    bulk_1 = FresnellSlab(n_1(λ), 2*π/λ, d_1, θ1)
    i_1_2_te, i_1_2_tm, θ2 = FresnellBoundrary(n_1(λ), n_2(λ), θ1)
    bulk_2 = FresnellSlab(n_2(λ), 2*π/λ, d_2, θ2)
    i_2_3_te, i_2_3_tm, θ3 = FresnellBoundrary(n_2(λ), n_3(λ), θ2)
    bulk_3 = FresnellSlab(n_3(λ), 2*π/λ, d_3, θ3)
    i_3_4_te, i_3_4_tm, θ4 = FresnellBoundrary(n_3(λ), n_4(λ), θ3)
    bulk_4 = FresnellSlab(n_4(λ), 2*π/λ, d_4, θ4)
    i_4_5_te, i_4_5_tm, θ5 = FresnellBoundrary(n_4(λ), n_5(λ), θ4)
    bulk_5 = FresnellSlab(n_5(λ), 2*π/λ, d_5, θ5)

    system_te = CascadeScattering([bulk_1, i_1_2_te, bulk_2, i_2_3_te, bulk_3, i_3_4_te, bulk_4, i_4_5_te, bulk_5])
    system_tm = CascadeScattering([bulk_1, i_1_2_tm, bulk_2, i_2_3_tm, bulk_3, i_3_4_tm, bulk_4, i_4_5_tm, bulk_5])

    te_r = abs(system_te[2, 1])^2
    te_t = abs(system_te[1, 1])^2
    tm_r = abs(system_tm[2, 1])^2
    tm_t = abs(system_tm[1, 1])^2

    [te_r, tm_r, te_t, tm_t][select]
end

function inside_slabs(expression)
    bulk_1(d) = FresnellSlab(n_1(λ), 2*π/λ, d, θ1)
    
    i_1_2_te, i_1_2_tm, θ2 = FresnellBoundrary(n_1(λ), n_2(λ), θ1)
    i_1_2_te_f(d) = i_1_2_te
    i_1_2_tm_f(d) = i_1_2_tm

    bulk_2(d) = FresnellSlab(n_2(λ), 2*π/λ, d, θ2)

    i_2_3_te, i_2_3_tm, θ3 = FresnellBoundrary(n_2(λ), n_3(λ), θ2)
    i_2_3_te_f(d) = i_2_3_te
    i_2_3_tm_f(d) = i_2_3_tm

    bulk_3(d) = FresnellSlab(n_3(λ), 2*π/λ, d, θ3)

    i_3_4_te, i_3_4_tm, θ4 = FresnellBoundrary(n_3(λ), n_4(λ), θ3)
    i_3_4_te_f(d) = i_3_4_te
    i_3_4_tm_f(d) = i_3_4_tm

    bulk_4(d) = FresnellSlab(n_4(λ), 2*π/λ, d, θ4)

    i_4_5_te, i_4_5_tm, θ5 = FresnellBoundrary(n_4(λ), n_5(λ), θ4)
    i_4_5_te_f(d) = i_4_5_te
    i_4_5_tm_f(d) = i_4_5_tm

    bulk_5(d) = FresnellSlab(n_5(λ), 2*π/λ, d, θ5)

    u_te, d_te = Interrogator(
        [bulk_1, i_1_2_te_f, bulk_2, i_2_3_te_f, bulk_3, i_3_4_te_f, bulk_4, i_4_5_te_f, bulk_5],
        [d_1, 0, d_2, 0, d_3, 0, d_4, 0, d_5],
        step,
        expression
    )
    u_tm, d_tm = Interrogator(
        [bulk_1, i_1_2_tm_f, bulk_2, i_2_3_tm_f, bulk_3, i_3_4_tm_f, bulk_4, i_4_5_tm_f, bulk_5],
        [d_1, 0, d_2, 0, d_3, 0, d_4, 0, d_5],
        step,
        expression
    )

    u_bp_te, u_bn_te, d_b_te = BorderInterrogator(
        [bulk_1(d_1), i_1_2_te, bulk_2(d_2), i_2_3_te, bulk_3(d_3), i_3_4_te, bulk_4(d_4), i_4_5_te, bulk_5(d_5)],
        [d_1, 0, d_2, 0, d_3, 0, d_4, 0, d_5],
        expression_borders
    )
    u_bp_tm, u_bn_tm, d_b_tm = BorderInterrogator(
        [bulk_1(d_1), i_1_2_tm, bulk_2(d_2), i_2_3_tm, bulk_3(d_3), i_3_4_tm, bulk_4(d_4), i_4_5_tm, bulk_5(d_5)],
        [d_1, 0, d_2, 0, d_3, 0, d_4, 0, d_5],
        expression_borders
    )

    u_te, d_te, u_tm, d_tm, u_bp_te, u_bn_te, d_b_te, u_bp_tm, u_bn_tm, d_b_tm
end


total_u_te, total_d_te, total_u_tm, total_d_tm, u_bp_te, u_bn_te, d_b_te, u_bp_tm, u_bn_tm, d_b_tm = inside_slabs(expression_total)
forward_u_te, forward_d_te, forward_u_tm, forward_d_tm, _, _, _, _, _, _ = inside_slabs(expression_forward)
backward_u_te, backward_d_te, backward_u_tm, backward_d_tm, _, _, _, _, _, _ = inside_slabs(expression_backward)

fig_inside = plot(title="Field inside " * description * "\nλ = $λ, θ = $θd", legend=:topleft, xlabel="distance in nm", ylabel="power", ticks=:native)
plot!(fig_inside, total_d_te .* 1e9, total_u_te, label="total te")
plot!(fig_inside, total_d_tm .* 1e9, total_u_tm, label="total tm")
plot!(fig_inside, forward_d_te .* 1e9, forward_u_te, label="forward te")
plot!(fig_inside, d_b_te .* 1e9, u_bp_te, seriestype =:scatter, label="forward te borders")
plot!(fig_inside, forward_d_tm .* 1e9, forward_u_tm, label="forward tm")
plot!(fig_inside, d_b_tm .* 1e9, u_bp_tm, seriestype =:scatter, label="forward tm borders")
plot!(fig_inside, backward_d_te .* 1e9, backward_u_te, label="backward te")
plot!(fig_inside, d_b_te .* 1e9, u_bn_te, seriestype =:scatter, label="backward te borders")
plot!(fig_inside, backward_d_tm .* 1e9, backward_u_tm, label="backward tm")
plot!(fig_inside, d_b_tm .* 1e9, u_bn_tm, seriestype =:scatter, label="backward tm borders")
display(fig_inside)
if save savefig(fig_inside, fig_path * "inside_" * fig_token * ".pdf") end


fig_angles = plot(title="By angle " * description * "\nλ = $λ", legend=:left, xlabel="angle in degrees", ylabel="power", ticks=:native)
plot!(fig_angles, θs .* (180/π), system_slabs.(θs, λ, 1), label="te reflection")
plot!(fig_angles, θs .* (180/π), system_slabs.(θs, λ, 2), label="tm reflection")
# plot!(fig_angles, θs .* (180/π), system_slabs.(θs, λ, 3), label="te transmission")
# plot!(fig_angles, θs .* (180/π), system_slabs.(θs, λ, 4), label="tm transmission")
display(fig_angles)
if save savefig(fig_angles, fig_path * "angles_" * fig_token * ".pdf") end

fig_wavelengths = plot(title="By wavelength " * description * "\nθ = $θd", legend=:right, xlabel="wavelength in nanometer", ylabel="power", ticks=:native)
plot!(fig_wavelengths, λs .* 1e9, system_slabs.(θ1, λs, 1), label="te reflection")
plot!(fig_wavelengths, λs .* 1e9, system_slabs.(θ1, λs, 2), label="tm reflection")
# plot!(fig_wavelengths, λs .* 1e9, system_slabs.(θ1, λs, 3), label="te transmission")
# plot!(fig_wavelengths, λs .* 1e9, system_slabs.(θ1, λs, 4), label="tm transmission")
display(fig_wavelengths)
if save savefig(fig_wavelengths, fig_path * "wavelengths_" * fig_token * ".pdf") end
