module fresnelltools

using LinearAlgebra

using ..matrixcore

export FresnellBoundrary, FresnellSlab


function FresnellBoundrary(n1::Real, n2::Real)::Matrix{Real}
    # Equation from Saleh & Teich 3.ed. eq.7.1-13
    # Returns the scattering matrix for the relevant border
    (1/(n1 + n2)) .* [
        2*n1  n2-n1;
        n1-n2 2*n2
    ]
end

function FresnellBoundrary(n1::Real, n2::Real, θ1)::(Matrix{Number}, Matrix{Number}, Number)
    # Equation from Saleh & Teich 3.ed. eq.7.1-23
    # It is similar to the other FresnellBoundrary function, but takes incident angle into account
    # Returns the TE scattering matrix, the TM scattering matrix and the outgoing transmitted angle

    # Handling total internal reflection
    if n1 > n2 & sin(θ1) * n1 / n2 > 1
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

    (Ste, Stm, θ2)
end

function FresnellSlab(n::Number, k_0::Number, d::Number)::Matrix{Number}
    # Equation from Saleh & Teich 3.ed. eq.7.1-4
    delay = ℯ^(-1im * n * k_0 * d)
    [
        delay 0;
        0im     delay
    ]
    
end


end # fresnelltools