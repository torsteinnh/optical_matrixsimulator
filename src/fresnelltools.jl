module fresnelltools

using LinearAlgebra

using ..matrixcore

export FresnellBoundrary, FresnellSlab, Grating, ContinousBorder, ThreeLayerSystem, Interrogator


function FresnellBoundrary(n1::Number, n2::Number)::Matrix{Number}
    # Equation from Saleh & Teich 3.ed. eq.6.2-8 & eq. 6.2-9
    # Returns the scattering matrix for the relevant border
    
    FresnellBoundrary(n1, n2, 0)[1]
end

function FresnellBoundrary(n1::Number, n2::Number, θ1::Number)::Tuple{Matrix{Number}, Matrix{Number}, Number}
    # Equation from Saleh & Teich 3.ed. eq.6.2-8 & eq. 6.2-9
    # It is similar to the other FresnellBoundrary function, but takes incident angle into account
    # Returns the TE scattering matrix, the TM scattering matrix and the outgoing transmitted angle

    # The imaginary coefficient 1e-14im has no physical reason to exist, and is only added for stability in the program.
    # It is assumed to be needed because of type instability in the asin function i Julia.
    # See https://github.com/JuliaLang/julia/issues/24296
    if ((typeof(n1) <: Complex) & (typeof(n2) <: Real))
        complex_modifier = 1e-14im
    else
        complex_modifier = 0im
    end
    θ2 = asin(complex_modifier + sin(θ1) * n1 / n2)

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
    interface2_te, interface2_tm, _ = FresnellBoundrary(n_bulk(λ), n_3, θ2)
    te_system = CascadeScattering([interface1_te, bulk, interface2_te])
    tm_system = CascadeScattering([interface1_tm, bulk, interface2_tm])

    te_system, tm_system
end

function Interrogator(layers::Vector{Function}, distances::Vector{Float64}, step::Float64)::Tuple{Vector{Float64}, Vector{Number}}
    # A small utility providing a default lambda for the Interrogator tool
    Power(Up, Um) = abs(Up + Um)^2
    Interrogator(layers, distances, step, Power)
end

function Interrogator(layers::Vector{Function}, distances::Vector{Float64}, step::Float64, expression::Function)::Tuple{Vector{Float64}, Vector{Number}}
    # A tool for inspecting the field inside a multilayer system

    total_system = Array{Matrix{Number}}(undef, length(layers))
    for i in 1:1:length(layers)
        total_system[i] = layers[i](distances[i])
    end
    U0_plus = 1
    U0_minus = (CascadeScattering(total_system) * [U0_plus, 0])[2]

    maxlen = ceil(Int, length(distances) * 2 + sum(distances) / step)
    d_values = Vector{Number}(undef, maxlen)
    u_values = Vector{Number}(undef, maxlen)

    d_values[1] = 0.0
    u_values[1] = expression(U0_plus, U0_minus)

    m_accumulated = I
    m_local = I
    m_partial = I
    d_accumulated = 0.0
    i = 1
    for (layer, distance) in zip(layers, distances)
        m_local = m_accumulated

        m_partial = StoM(layer(step))
        for _ in  0:step:distance
            i += 1
            d_accumulated += step
            m_local = m_partial * m_local

            U_plus, U_minus = m_local * [U0_plus, U0_minus]
            measure = expression(U_plus, U_minus)
            d_values[i] = d_accumulated
            u_values[i] = measure
        end
        
        d_rest = distance % step
        if d_rest != 0
            m_partial = StoM(layer(d_rest))
            i += 1
            d_accumulated += d_rest
            m_local = m_partial * m_local

            (U_plus, U_minus) = m_local * [U0_plus, U0_minus]
            
            measure = expression(U_plus, U_minus)
            d_values[i] = d_accumulated
            u_values[i] = measure
        end

        m_accumulated = StoM(layer(distance)) * m_accumulated
    end

    u_values[1:i], d_values[1:i]
end


end # fresnelltools