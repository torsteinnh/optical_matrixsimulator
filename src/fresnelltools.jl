module fresnelltools

using LinearAlgebra

using ..matrixcore

export FresnellBoundrary, FresnellSlab, Grating, ContinousBorder, ThreeLayerSystem


function FresnellBoundrary(n1::Number, n2::Number)::Matrix{Number}
    # Equation from Saleh & Teich 3.ed. eq.6.2-8 & eq. 6.2-9
    # Returns the scattering matrix for the relevant border
    
    FresnellBoundrary(n1, n2, 0)[1]
end

function FresnellBoundrary(n1::Number, n2::Number, θ1::Number)::Tuple{Matrix{Number}, Matrix{Number}, Number}
    # Equation from Saleh & Teich 3.ed. eq.6.2-8 & eq. 6.2-9
    # It is similar to the other FresnellBoundrary function, but takes incident angle into account
    # Returns the TE scattering matrix, the TM scattering matrix and the outgoing transmitted angle

    θ2 = asin(0im + sin(θ1) * n1 / n2)

    r_te_12 = (n1 * cos(θ1) - n2 * cos(θ2)) / (n1 * cos(θ1) + n2 * cos(θ2))
    r_te_21 = (n2 * cos(θ2) - n1 * cos(θ1)) / (n2 * cos(θ2) + n1 * cos(θ1))

    r_tm_12 = (n1 * sec(θ1) - n2 * sec(θ2)) / (n1 * sec(θ1) + n2 * sec(θ2))
    r_tm_21 = (n2 * sec(θ2) - n1 * sec(θ1)) / (n2 * sec(θ2) + n1 * sec(θ1))

    t_te_12 = 1 + r_te_12
    t_te_21 = 1 + r_te_21

    t_tm_12 = (1 + r_tm_12) * cos(θ1) / cos(θ2)
    t_tm_21 = (1 + r_tm_21) * cos(θ2) / cos(θ1)


    Ste = [
        t_te_12 r_te_21;
        r_te_12 t_te_21
    ]
    Stm = [
        t_tm_12 r_tm_21;
        r_tm_12 t_tm_21
    ]

    return Ste, Stm, θ2
end

function FresnellSlab(n::Number, k_0::Number, d::Number, θ::Number)::Matrix{Number}
    # Equation from Saleh & Teich 3.ed. eq.7.1-4
    delay = ℯ^(-1im * n * k_0 * d * cos(θ))
    S = [
        delay 0;
        0     delay
    ]

    if abs(S[1, 1]) < 1e-5
        S /= abs(S[1, 1])
        S *= 1e-5
    end

    S
end

function Grating(n_spechial::Number, n_normal::Number, d_spechial::Real, d_normal::Real, layers::Int, k_0::Number)::Matrix{Number}
    # A utility for creating the scattering matrix of a periodic grating

    border_normal_spechial = FresnellBoundrary(n_normal, n_spechial)
    bulk_spechial = FresnellSlab(n_spechial, k_0, d_spechial, 0)
    border_spechial_normal = FresnellBoundrary(n_spechial, n_normal)
    bulk_normal = FresnellSlab(n_normal, k_0, d_normal, 0)

    cell = CascadeScattering([border_normal_spechial, bulk_spechial, border_spechial_normal, bulk_normal])
    Repeating(cell, layers)
end

function Grating(n_spechial::Number, n_normal::Number, d_spechial::Real, d_normal::Real, layers::Int, k_0::Number, θ1::Real)::Tuple{Matrix{Number}, Matrix{Number}}
    # A utility similar to Grating, but it takes an angle and assumes an infinitely wide grating
    # Returns the TE and TM scattering matrices for the grating
    # Unlike the angleded FresnellBoundrary tool, no angle is returned, this is because the grating assumes normal material on both sides

    n_s_te, n_s_tm, θ2 = FresnellBoundrary(n_normal, n_spechial, θ1)
    s_n_te, s_n_tm, _  = FresnellBoundrary(n_spechial, n_normal, θ2)

    n_b = FresnellSlab(n_normal, k_0, d_normal, θ1)
    s_b = FresnellSlab(n_spechial, k_0, d_spechial, θ2)

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
        components[2 * i] = FresnellSlab(n2, k_0, d_stepp, 0)
    end

    CascadeScattering(components)
end

function ThreeLayerSystem(n_1::Number, n_bulk, n_3::Number, λ, θ, d)
    interface1_te, interface1_tm, θ2 = FresnellBoundrary(n_1, n_bulk(λ), θ)
    bulk = FresnellSlab(n_bulk(λ), 2*π / λ, d, θ2)
    interface2_te, interface2_tm, _ = FresnellBoundrary(n_bulk(λ), n_3, θ2 + 1e-14im)
    te_system = CascadeScattering([interface1_te, bulk, interface2_te])
    tm_system = CascadeScattering([interface1_tm, bulk, interface2_tm])

    te_system, tm_system
end


end # fresnelltools