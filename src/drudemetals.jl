module drudemetals

using ..utilities

export θ_plasmone, n_drude


function ϵ(n::Number)::Number
    real(n)^2 - imag(n)^2 + 1im * 2 * real(n) * imag(n)
end

function n_drude(σ_0, τ, λ)
    ω = 2 * π * c_0 / λ
    √ϵ_c_drude(σ_0, τ, ω)
end


function ϵ_c_drude(σ_0, τ, ω)
    # jsee jsaleh & Teich 3.ed eq.8.2-11

    ω_p = √(σ_0 / (ϵ_0 * τ))
    ϵ_0 * (1 + ω_p^2 / (-(ω^2) + 1im * ω / τ))
end


function β_plasmone(n1, n2, k_0)
    # See Maier eq. 2.14 and eq. 1.11a and b

    ϵ_1 = ϵ(n1)
    ϵ_2 = ϵ(n2)
    k_0 * √(ϵ_1 * ϵ_2 / (ϵ_1 + ϵ_2))
end

function θ_plasmone(ni, n1, n2, k_0)
    β = real(β_plasmone(n1, n2, k_0))
    ki = real(ni) * k_0

    if ((β / ki) > 1)
        return NaN
    end
    asin(β / ki)
end


end # drudemetals
