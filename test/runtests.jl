using AnimatedOptimization
using Test


@testset "Escape local minima" begin
    # Test heuristic optimization with function with many local minima, one global (creased egg carton)
    using LinearAlgebra, Distributions

    function creased_egg_carton(A, ζ, w0, φ)
      @assert ζ <= 1.0
      t -> (A*(1-exp(-ζ*w0*abs(t[1]-2)))+sin(sqrt(1-ζ^2)*w0*abs(t[1]-2)+φ) + A*(1-exp(-ζ*w0*abs(t[2]-1)))+sin(sqrt(1-ζ^2)*w0*abs(t[2]-1)+φ) )^1
    end
    f = creased_egg_carton(10, 0.5, 5, 0.0)

    # Solver parameters
    x0 = [-0.5, -0.5]
    vtol = 1E-3

    # Solve and test
    sol_mrs = minrandomsearch(f, x0, 20, var0=1, vshrink=0.5, vtol=vtol )
    @test norm(sol_mrs[2] - [2.0 , 1.0], Inf) < vtol

end
