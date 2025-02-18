using Test
using DiffMinPoly
using StructuralIdentifiability

include("test_utils.jl")

cases = []

generics = [
    # [2, 2],
    # [1, 2, 2],
    [2, 1, 1],
    [2, 2, 2],
    [2, 3, 3],
    # [2, 4, 4], 
    # [2, 5, 5],  
    # [3, 1, 1], 
    # [3, 2, 2],
    # [3, 3, 3],
    # [1,2,2,2],
    # [2,1,1,1],  
]

for ds in generics
    println(ds)
    push!(cases, rand_ode(ds))
end

push!(
    cases,
    @ODEmodel(
        x1'(t) = x1(t) + 3 * x1(t) * x2(t),
        x2'(t) = -5 * x2(t) + x1(t) * x2(t),
        y(t) = x1(t)
    ),

    @ODEmodel(
       m'(t) = 4*m(t) + 3 * n(t)^3 * m(t) + n(t),               
       n'(t) = -2 * n(t) + -5*m(t) * n(t) + m(t)^2,                
       y(t) = m(t)
    ),


    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t),
        x2'(t) = -3* x2(t) + 2*x1(t) * x2(t),
        y(t) = x1(t)
    ),

    
    @ODEmodel(
       x1'(t) = 2 * x1(t)^2 + 17*x2(t)^2 + 3*x3(t)^2 + 1,
       x2'(t) = x2(t)^2,
       x3'(t) = x3(t)^2,
       y(t) = x1(t)
    ),

    @ODEmodel(
       x1'(t) = 2^5 * x1(t)^2 - 3^4 * x2(t) + 7 - 9 * x3(t)*x2(t),             
       x2'(t) = 2^3 * x1(t) + x3(t)^2 + 5 + x2(t)^2,             
       x3'(t) = -2 * x3(t)*x2(t) + 2^4 * x1(t) * x2(t)^2,             
       y(t) = x1(t)
    ),
    
    @ODEmodel(
       x1'(t) = 2^31 * x1(t)^2 + 2^32 *x2(t)^2 + 3^18 *x3(t)^2 + 2^30,
       x2'(t) = 2^31 * x2(t)^2 + 2^15 * x1(t),
       x3'(t) = 2^12 * x3(t)^2 + 2^12 * x1(t) + 2^11 * x2(t)^2,
       y(t) = x1(t)
    ),
    
    # @ODEmodel(
    #    x1'(t) = 3 * x1(t)^2 + 16 * x2(t)^2 + 18 * x1(t) + 42,
    #    x2'(t) = 2^31 * x2(t)^2 + 2^15 * x1(t) + 17,
    #    x3'(t) = 19 * x1(t) + 95 * x2(t),
    #    y(t) = x1(t)
    # ),
)

@testset "Testing against the standard algorithms" begin
    for c in cases
        @info c
        @test check_ansatz(c)
    end
end
