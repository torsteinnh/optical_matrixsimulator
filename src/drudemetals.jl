module drudemetals

using ..utilities

export σ_drude, n_metal, θ_plasmone


function ϵ(n::Number)::Number
    real(n)^2 - imag(n)^2 + 1im * 2 * real(n) * imag(n)
end


function σ_drude(σ_0, ω, τ)
    # See Saleh & Teich 3.ed eq.8.2-10

    σ_0 / (1 + 1im * ω * τ)
end

function n_metal(σ, γ, ϵ, ω)::Number
    # See Saleh & Teich 3.ed eq.8.2-6

    √(ϵ / ϵ_0) * √(1 + σ / (1im * ω * ϵ)) + 1im * γ * c_0 / ω
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
