module analyticalmaterials

using ..utilities

export θ_plasmone, n_drude, n_MaxwellGarnett


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
    # jsee jsaleh & Teich 3.ed eq.8.2-11

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


function n_MaxwellGarnett(ϵ_a::Number, ϵ_b::Number, f_A::Real)::Number
    # Hans Arwin, eq. 2.23
    # Coated spheres model
    # f_A is the fill factor

    ϵ_total = ϵ_b * (ϵ_a + 2*ϵ_b + 2*f_A*(ϵ_a - ϵ_b)) / (ϵ_a + 2*ϵ_b - f_A*(ϵ_a - ϵ_b))
    n(ϵ_total)
end


end # analyticalmaterials
