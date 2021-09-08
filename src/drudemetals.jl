module drudemetals

export σ_drude


function σ_drude(σ_0, ω, τ)
    # See Saleh & Teich 3.ed eq.8.2-10

    σ_0 / (1 + 1im * ω * τ)
end


end # drudemetals
