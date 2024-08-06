include("solver_love_and_support.jl")
include("test_utils.jl")

using Test

cases = []

generics = [
    [1, 1, 1],
    [1, 2, 3],
    [1, 2, 2],
    [1, 1, 2],
    [2, 1, 1],
    [2, 2, 2],
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

@testset "Testing against the standard algorithms" begin
    for c in cases
        @test check_ansatz(c)
    end
end
