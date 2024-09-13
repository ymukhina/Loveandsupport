include("solver_love_and_support.jl")
include("test_utils.jl")

using Test

cases = []

generics = [
    [2, 1, 1],
    [2, 2, 2],
    [2, 3, 3],
    [2, 4, 4], 
    [2, 5, 5],  
    [3, 1, 1], 
    [3, 2, 2],
    [3, 3, 3],
    [1,2,2,2],
    [2,1,1,1],  
]

for ds in generics
    push!(cases, rand_ode(ds))
end

push!(
    cases,
    @ODEmodel(
        x1'(t) = x1(t) + 3 * x1(t) * x2(t),
        x2'(t) = -5 * x2(t) + x1(t) * x2(t),
        y(t) = x1(t)
    )
)

# Gleb: the following test cases should be added:
#   * where the result of elimination is of lower order
#   * where the coefficients are very large
#   * the one which was a counter example for the power series approach

@testset "Testing against the standard algorithms" begin
    for c in cases
        @test check_ansatz(c)
    end
end