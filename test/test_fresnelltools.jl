using Test

include("../src/simulator.jl")
using .simulator.fresnelltools


@testset "test_fresnelltools" begin
    
    @testset "FresnellBoundraryNoangle" begin
      # A randomized test verifying the angle-independent fresnell boundrary
      # is equivalent to the angle dependent at 90 degrees.
      for _ in 1:100
        n1 = rand(ComplexF64)
        n2 = rand(ComplexF64)

        independent = FresnellBoundrary(n1, n2)
        te_matrix, tm_matrix, θ = FresnellBoundrary(n1, n2, 0)

        @test independent == te_matrix == tm_matrix
        @test θ == 0
      end
    end

    @testset "FresnellBoundraryNonexisting" begin
      # Testing a variety of angles and refractive indices to check that nonexisting borders do no reflect.
      for _ in 1:100
        θ1 = rand() * π / 2
        n = rand() * 2
        te, tm, θ2 = FresnellBoundrary(n, n, θ1)

        @test θ1 ≈ θ2
        @test round.(te, digits=5) ≈ round.(tm, digits=5) ≈ [1.0 0.0; 0.0 1.0]
      end
    end

    @testset "FresnellSlabTest" begin
      for _ in 1:100
        length = rand() * 10
        slab = FresnellSlab(1, 1, length, 0)

        @test slab ≈ ℯ^(- length * 1e-6) .* [ℯ^(- length * 1im) 0; 0 ℯ^(- length * 1im)]
      end

      @test FresnellSlab(1, 2, 0, 0) == [1 0; 0 1]
    end

end
