module matrixcore

export StoM, MtoS, CascadeScattering, Repeating

using LinearAlgebra


function StoM(S::Matrix{T})::Matrix{<:Number} where T <: Number
    # Converts a 2x2 scattering matrix to a transfer matrix
    # See Saleh & Teich 3.ed eq.7.1-6
    # Verified manually
    @assert(size(S) == (2, 2))
    @inbounds begin
        t12 = S[1, 1]
        r21 = S[1, 2]
        r12 = S[2, 1]
        t21 = S[2, 2]
    end
        M = [
        t12*t21 - r12*r21 r21;
        -r12              1
    ]
    (1/t21) .* M
end

function MtoS(M::Matrix{T})::Matrix{<:Number} where T <: Number
    # Converts a 2x2 transfer matrix to a scattering matrix
    # See Saleh & Teich 3.ed eq.7.1-5
    # Verified by tests as inverfse of StoM
    StoM(M)
end


function CascadeScattering(layers::Vector)::Matrix{Number}
    # Generates the total S matrix for a layered system
    # Matrixes are given in the order the elements appear, the array is inversed inside the function
    accumulated = I

    for S in view(layers, length(layers):-1:1)
        M = StoM(S)
        accumulated *= M
    end

    MtoS(accumulated)
end


function Repeating(cell::Matrix{T}, repetition::Number)::Matrix{T} where T <: Number
    # Simple helpertool to make sure gratings are handleded correctly
    # This helper helps achieve the goal that transfer matrices are only handled in this module
    M = StoM(cell)
    MtoS(M ^ repetition)
end


end # matrixcore
