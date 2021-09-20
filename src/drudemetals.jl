module drudemetals

using ..utilities

export σ_drude, n_metal


function σ_drude(σ_0, ω, τ)
    # See Saleh & Teich 3.ed eq.8.2-10

    σ_0 / (1 + 1im * ω * τ)
end

function n_metal(σ, γ, ϵ, ω)::Number
    # See Saleh & Teich 3.ed eq.8.2-6

    √(ϵ / ϵ_0) * √(1 + σ / (1im * ω * ϵ)) + 1im * γ * c_0 / ω
end


end # drudemetals
