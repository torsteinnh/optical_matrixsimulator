module analyticalmaterials

using ..utilities

export θ_plasmone, n_drude, θ_grating_coupling, ϵ, n, hs_Palm


function ϵ(n::Number)::Number
    real(n)^2 - imag(n)^2 + 1im * 2 * real(n) * imag(n)
end

function n(ϵ::Number)::Number
    √ϵ
end

function n_drude(σ_0::Real, τ::Real, λ::Real)::Number
    ω = 2 * π * c_0 / λ
    √ϵ_c_drude(σ_0, τ, ω)
end


function ϵ_c_drude(σ_0::Real, τ::Real, ω::Real)::Number
    # See Saleh & Teich 3.ed eq.8.2-11

    ω_p = √(σ_0 / (ϵ_0 * τ))
    ϵ_0 * (1 + ω_p^2 / (-(ω^2) + 1im * ω / τ))
end


function β_plasmone(n1::Number, n2::Number, k_0::Number)::Number
    # See Maier eq. 2.14 and eq. 1.11a and b

    ϵ_1 = ϵ(n1)
    ϵ_2 = ϵ(n2)
    k_0 * √(ϵ_1 * ϵ_2 / (ϵ_1 + ϵ_2))
end

function θ_plasmone(ni::Number, n1::Number, n2::Number, k_0::Number)::Number
    β = real(β_plasmone(n1, n2, k_0))
    ki = real(ni) * k_0

    if ((β / ki) > 1)
        return NaN
    end
    asin(β / ki)
end


function θ_grating_coupling(λ, neff_core, n_cladding, Λ)
    asin(((λ * 2 / Λ) - neff_core) / n_cladding)
end

function θ_grating_coupling(λ, n_core)
    # A default for the grating in question
    θ_grating_coupling(λ, 1.4676, n_core, 1073.81e-9)
end


module hs_Palm
# An approximation of the nonlinear function h(c) such that ϵ(c) = h(c)ϵ(0) for palladium and gold-palladium alloys.
# Estimate based on linear approximation from Palm2019 and theory form Tobiska2001.
# This module should only be used as a first approximation, not for spesific optimization.

using Interpolations
using DelimitedFiles
using Match

import ..ϵ


export h


rawdatafile034 = readdlm("materials/palm_PdAu-H/HM_Pd034.csv", ',', skipstart=2)
rawdatafile042 = readdlm("materials/palm_PdAu-H/HM_Pd042.csv", ',', skipstart=2)
rawdatafile052 = readdlm("materials/palm_PdAu-H/HM_Pd052.csv", ',', skipstart=2)
rawdatafile073 = readdlm("materials/palm_PdAu-H/HM_Pd073.csv", ',', skipstart=2)
rawdatafile100 = readdlm("materials/palm_PdAu-H/HM_Pd100.csv", ',', skipstart=2)


function loadLine(filerawdata, λ::Float64, Pd_c::Float64)::Function
    point = @match λ begin
        300e-9 => 1
        600e-9 => 2
        900e-9 => 3
        1200e-9 => 4
        1500e-9 => 5
        _ => 0
    end
    offsett = 10

    raw_n = filerawdata[:, (point*2 - 1):(point*2)]
    raw_k = filerawdata[:, (point*2 - 1 + offsett):(point*2 + offsett)]

    raw_n = raw_n[(raw_n[:, 1] .!= ""), :]
    raw_k = raw_k[(raw_k[:, 1] .!= ""), :]

    sorted_n = raw_n[sortperm(raw_n[:, 1]), :]
    sorted_k = raw_k[sortperm(raw_k[:, 1]), :]

    n_estimator = LinearInterpolation(sorted_n[:, 1], sorted_n[:, 2])
    k_estimator = LinearInterpolation(sorted_k[:, 1], sorted_k[:, 2])

    h_total(c) = ϵ(n_estimator(c_to_HM(c, Pd_c)) + k_estimator(c_to_HM(c, Pd_c)) * 1im) / ϵ(sorted_n[1, 2] + sorted_k[1, 2] * 1im)

    h_total
end


function c_to_HM(H_c::Float64, Pd_c::Float64)::Float64
    # Based on Palm2019 fig. S7
    # Assumes concentration under 1 atmosphere pressure
    # Assumes the effect of a low pressure of pure hydrogen is the same as the equivalent partial pressure hydrogen in an otherwize inert atmosphere totalling 1 atmosphere
    # Assumes the relation between HM and H_c is proportional to the square root of H_c, as shown in separate note
    # H_c is the partial pressure of hydrogen in the atmosphere given in atmospheres, the approximation is only valid for H_c ⪅ 0.25
    # Pd_c is the atomic alloying ratio between Pd and Au, and should be between 0.4 and 1    
    # The following data was extracted from the figure:
    # 41.50974025974026	0.0011235955056179137
    # 99.9512987012987	0.6089887640449437
    
    HM_Pd_c025 = 1.04012484394507 * Pd_c - 0.43062952559301

    sqrt(H_c * 1.01325) * HM_Pd_c025 / sqrt(0.25)
end

function h(λ::Float64, H_c::Float64, Pd_c::Float64)::Number
    # H_c and Pd_c is the hydrogen concentration and the Pd atomic concentration in the Pd-Au alloy, both between 0 and 1
    
    filerawdata = @match Pd_c begin
        0.34 => rawdatafile034
        0.42 => rawdatafile042
        0.52 => rawdatafile052
        0.73 => rawdatafile073
        1.0 => rawdatafile100
        _ => "Illegal Pd_c"
    end
    
    λ_lower = 300e-9
    λ_upper = 1500e-9
    if 300e-9 <= λ <= 600e-9
        λ_lower = 300e-9
        λ_upper = 600e-9
    elseif 600e-9 <= λ <= 900e-9
        λ_lower = 600e-9
        λ_upper = 900e-9
    elseif 900e-9 <= λ <= 1200e-9
        λ_lower = 900e-9
        λ_upper = 1200e-9
    elseif 1200e-9 <= λ <= 1500e-9
        λ_lower = 1200e-9
        λ_upper = 1500e-9
    end

    h_lower = loadLine(filerawdata, λ_lower, Pd_c)(H_c)
    h_upper = loadLine(filerawdata, λ_upper, Pd_c)(H_c)

    λ_partial = (λ - λ_lower) / (λ_upper - λ_lower)
    h_real = (h_upper - h_lower) * λ_partial + h_lower

    h_real
end

end


end # analyticalmaterials
