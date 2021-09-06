using Test

include("../src/matrixcore.jl")
using .matrixcore


@testset "test_matrixcore" begin
    
    @testset "StoMtoS" begin
        # Running random tests to check that MtoS and StoM are each others inverse
        for _ in 1:100
            S = rand(ComplexF64, 2, 2)
            @test S ≈ (StoM ∘ MtoS)(S)
        end
    end

    @testset "test_cascades" begin
        # Getting random transmittance and reflectances, testing two mirror systems
        # Saleh & Teich 3.ed eq.7.1-8
        for _ in 1:100
            t12 = rand(ComplexF64)
            t21 = rand(ComplexF64)

            r12 = rand(ComplexF64)
            r21 = rand(ComplexF64)

            t23 = rand(ComplexF64)
            t32 = rand(ComplexF64)

            r23 = rand(ComplexF64)
            r32 = rand(ComplexF64)

            ϕ = rand() * 2π

            mirror1 = [
                t12 r21;
                r12 t21
            ]
            air = [ℯ^(-ϕ * im) 0; 0 ℯ^(-ϕ * im)]
            mirror2 = [
                t23 r32;
                r23 t32
            ]

            machine = CascadeScattering([mirror1, air, mirror2])

            @test machine[1, 1] ≈ (t12 * t23 * ℯ^(-ϕ * im)) / (1 - (r21 * r23 * ℯ^(-ϕ * 2im)))
            @test machine[2, 1] ≈ r12 + (t12 * t21 * r23 * ℯ^(-ϕ * 2im)) / (1 - (r21 * r23 * ℯ^(-ϕ * 2im)))
        end
    end
    
end
