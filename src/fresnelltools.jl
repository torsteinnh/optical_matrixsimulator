module fresnelltools

using LinearAlgebra

using ..matrixcore

export FresnellBoundrary, FresnellSlab, Grating, ContinousBorder, PlasmoneScattering


function FresnellBoundrary(n1::Number, n2::Number)::Matrix{Number}
    # Equation from Saleh & Teich 3.ed. eq.7.1-13
    # Returns the scattering matrix for the relevant border
    (1/(n1 + n2)) .* [
        2*n1  n2-n1;
        n1-n2 2*n2
    ]
end

function FresnellBoundrary(n1::Number, n2::Number, θ1::Number)::Tuple{Matrix{Number}, Matrix{Number}, Number}
    # Equation from Saleh & Teich 3.ed. eq.7.1-23
    # It is similar to the other FresnellBoundrary function, but takes incident angle into account
    # Returns the TE scattering matrix, the TM scattering matrix and the outgoing transmitted angle

    # Handling total internal reflection
    θ2 = 0
    if ((abs(n1) > abs(n2)) && (abs((sin(θ1) * n1 / n2)) > 1))
        θ2 = asin(0im + sin(θ1) * n1 / n2)
    else
        θ2 = asin(sin(θ1) * n1 / n2)
    end

    ten1 = n1 * cos(θ1)
    ten2 = n2 * cos(θ2)

    tmn1 = n1 * sec(θ1)
    tmn2 = n2 * sec(θ2)
    tma = cos(θ1) / cos(θ2)

    Ste = FresnellBoundrary(ten1, ten2)
    Stm = FresnellBoundrary(tmn1, tmn2)
    Stm[1, 1] *= tma
    Stm[2, 2] /= tma

    return Ste, Stm, θ2
end

function FresnellSlab(n::Number, k_0::Number, d::Number)::Matrix{Number}
    # Equation from Saleh & Teich 3.ed. eq.7.1-4
    delay = ℯ^(-1im * n * k_0 * d)
    dampening = ℯ^(- d * 1e-6)
    dampening .* [
        delay 0;
        0     delay
    ]
end

function Grating(n_spechial::Number, n_normal::Number, d_spechial::Real, d_normal::Real, layers::Int, k_0::Number)::Matrix{Number}
    # A utility for creating the scattering matrix of a periodic grating

    border_normal_spechial = FresnellBoundrary(n_normal, n_spechial)
    bulk_spechial = FresnellSlab(n_spechial, k_0, d_spechial)
    border_spechial_normal = FresnellBoundrary(n_spechial, n_normal)
    bulk_normal = FresnellSlab(n_normal, k_0, d_normal)

    cell = CascadeScattering([border_normal_spechial, bulk_spechial, border_spechial_normal, bulk_normal])
    Repeating(cell, layers)
end

function Grating(n_spechial::Number, n_normal::Number, d_spechial::Real, d_normal::Real, layers::Int, k_0::Number, θ1::Real)::Tuple{Matrix{Number}, Matrix{Number}}
    # A utility similar to Grating, but it takes an angle and assumes an infinitely wide grating
    # Returns the TE and TM scattering matrices for the grating
    # Unlike the angleded FresnellBoundrary tool, no angle is returned, this is because the grating assumes normal material on both sides

    n_s_te, n_s_tm, θ2 = FresnellBoundrary(n_normal, n_spechial, θ1)
    s_n_te, s_n_tm, _  = FresnellBoundrary(n_spechial, n_normal, θ2)

    n_b = FresnellSlab(n_normal, k_0, d_normal / cos(θ1))
    s_b = FresnellSlab(n_spechial, k_0, d_spechial / cos(θ2))

    cell_te = CascadeScattering([n_s_te, s_b, s_n_te, n_b])
    cell_tm = CascadeScattering([n_s_tm, s_b, s_n_tm, n_b])

    Repeating(cell_te, layers), Repeating(cell_tm, layers)
end

function ContinousBorder(n_from::Number, n_to::Number, d::Real, k_0::Number, stepps::Int)::Matrix{Number}
    # A utility for creating a semi-continous change in refractive index
    # The gradient is approximated as linear

    components = Array{Matrix{Number}}(undef, 2 * stepps)

    n_stepp = (n_to - n_from) / stepps
    d_stepp = d / stepps
    for i in 1:stepps
        n1 = n_from + (i - 1) * n_stepp
        n2 = n1 + n_stepp

        components[2 * i - 1] = FresnellBoundrary(n1, n2)
        components[2 * i] = FresnellSlab(n2, k_0, d_stepp)
    end

    CascadeScattering(components)
end

function PlasmoneScattering(β, β0, Δβ)::Number
    # A model for reflection coefficient R at a prism-plasmone-coupling
    # Equation from Maier(2007) eq. 3.1
    # β is dependent on the SPP system, β0 is the "Singel interface value"
    # β = β0 + Re{Δβ}

    Γlr = imag(Δβ)
    Γabs = imag(β0)

    1 - (4 * Γlr * Γabs) / ((β - (β0 + Δβ))^2 + (Γlk + Γabs)^2)
end


end # fresnelltools