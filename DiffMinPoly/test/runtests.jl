using Test
using DiffMinPoly
using StructuralIdentifiability

include("test_utils.jl")

cases = []

generics = [
    [1, 2], 
    [1, 1], 
    [2, 1],
    [2, 2],
    [1, 2, 2],
    [2, 1, 1],
    # [2, 2, 1],
    # [2, 2, 2],
    # [2,2,2,1],
    # [2, 3, 3],
    # [2, 4, 4], 
    # [2, 5, 5],  
    # [3, 1, 1], 
    # [3, 2, 2],
    # [3, 3, 3],
    # [1,2,2,2],
    # [2,1,1,1],  
]

for ds in generics
    push!(cases, rand_ode_x(ds))
    push!(cases, rand_ode(ds))
   # push!(cases, rand_ode([3,2,2,1]))
end

push!(
    cases,

    @ODEmodel(
        x1'(t) = x1(t) + 3 * x1(t) * x2(t),
        x2'(t) = -5 * x2(t) + x1(t) * x2(t),
        y(t) = x1(t)^2 + x2(t)
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t),
        x2'(t) = -3 * x2(t) + 2*x1(t) * x2(t),
        y(t) = x1(t)
    ),

    @ODEmodel(
        x1'(t) = x1(t)^2+ x2(t) + 1,
        x2'(t) = 2* x1(t) + 3 * x2(t) + 4,
        y(t) = x1(t) + x2(t) + 1
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t)^2,
        x2'(t) = -3 * x2(t) + 2*x1(t) * x2(t),
        y(t) = x2(t) + x1(t)
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t)^2,
        x2'(t) = -3 * x2(t)^2 + 2*x1(t) * x2(t),
        y(t) = x2(t) + x1(t)
    ),


    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t)^2,
        x2'(t) = -3 * x2(t)^2 + 2*x1(t) * x2(t),
        y(t) = 3*x2(t) + 45*x1(t) 
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t)^3,
        x2'(t) = -3 * x2(t)^3 + 2*x1(t) * x2(t),
        y(t) = 2*x2(t) + 5*x1(t) 
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t)^2,
        x2'(t) = -3 * x2(t)^2 + 2*x1(t) * x2(t),
        x3'(t) = 6 * x2(t)^2 + 9*x1(t) * 32*x2(t) + 16*x1(t),
        y(t) = x1(t) + 7*x2(t) 
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t)^3 - x2(t)^2,
        x2'(t) = -3 * x2(t)^2 + 2*x1(t) * x2(t),
        x3'(t) = 6 * x2(t)^2 + 9*x1(t) * 32*x2(t) + 16*x1(t),
        y(t) = x1(t) + 5 * x2(t) 
    ),
    

    @ODEmodel(
        x1'(t) = 3*x3(t) - x2(t),
        x2'(t) = -3 * x2(t) + 2*x1(t) * x2(t),
        x3'(t) = 6 * x1(t) +  32*x3(t),
        y(t) = x1(t) + 35 * x3(t) + x2(t)
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t)^3 - x2(t)^2,
        x2'(t) = -3 * x2(t)^2 + 2*x1(t) * x2(t),
        x3'(t) = 6 * x2(t)^2 + 9*x1(t) * 32*x2(t) + 16*x1(t),
        y(t) = x1(t) + 5 * x2(t) + 1
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t)^2,
        x2'(t) = -3 * x2(t)^2 + 2*x1(t) * x2(t),
        y(t) = 3*x2(t) + 3 
    ),


    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t)^2,
        x2'(t) = -3 * x2(t)^2 + 2*x1(t) * x2(t),
        y(t) = 3*x2(t) + 45*x1(t)  + 13
    ),

    @ODEmodel(
        x1'(t) = 3*x3(t)^2 - x2(t),
        x2'(t) = -3 * x2(t) + 2*x1(t) * x2(t),
        x3'(t) = 6 * x1(t)^2 +  32*x3(t),
        y(t) = x1(t) + 35 * x3(t) + 6
    ),

    @ODEmodel(
        x1'(t) = 3*x3(t)^2 - x2(t),
        x2'(t) = -3 * x2(t) + 2*x1(t),
        x3'(t) = 6 * x1(t) +  32*x3(t),
        y(t) = x1(t)
    ),

    @ODEmodel(
        x1'(t) = 3*x3(t)^2 - x2(t) + x1(t),
        x2'(t) = -3 * x2(t) + 2*x1(t),
        x3'(t) = 6 * x1(t) +  32*x3(t),
        y(t) = x1(t) + 35 * x3(t) + 6
    ),

    @ODEmodel(
        x1'(t) = 3*x3(t)^2 - x2(t) + x1(t),
        x2'(t) = -3 * x2(t) + 2*x1(t)^2,
        x3'(t) = 6 * x1(t) +  32*x3(t)^2,
        y(t) = x1(t) + 35 * x3(t) + 6
    ),


    
    @ODEmodel(
        x1'(t) = 3*x3(t)^2 - x2(t) + x1(t),
        x2'(t) = - 3 * x2(t) + 2*x1(t)^2,
        x3'(t) = 6 * x1(t) +  32*x3(t)^2,
        y(t) = x1(t) + 35 * x2(t) + 6 
    ),


    
     @ODEmodel(
        x1'(t) = 3*x3(t)^2 - x2(t) + x1(t) + 35 * (-3 * x2(t) + 2*x1(t)^2),
        x2'(t) = -3 * x2(t) + 2*(x1(t) - 35 * x2(t) - 6)^2,
        x3'(t) = 6 * (x1(t) - 35 * x2(t) - 6)  +  32*x3(t)^2,
        y(t) = x1(t) 
    ),


    @ODEmodel(
        x1'(t) =  -x2(t) + 3*x3(t)^2,
        x2'(t) =2*x1(t) - 3*x2(t),
        x3'(t) = 6*x1(t) + 32*x3(t),
        y(t) = x1(t)

    ),

    @ODEmodel(
        x1'(t) =  -x2(t) + 3*x1(t)^2,
        x2'(t) =2*x1(t)^2 - 3*x2(t),
        y(t) = x1(t)*x2(t) + x1(t) + 7 * x2(t)

    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) + 7*x2(t) - 6,
        x2'(t) = -9*x1(t)^2 + x1(t)*x2(t) - 5*x1(t) + 2*x2(t)^2 + 4*x2(t) - 8,
        y(t) = 3*x1(t)^2 - 4*x1(t)*x2(t) - 2*x1(t) - 7*x2(t)^2 - 6*x2(t) + 9
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) + 7*x2(t) - 6,
        x2'(t) = 9*x1(t) + 17*x2(t) - 16,
        y(t) = x1(t)^2 - (x2(t)^2 - x2(t)^3)
    ),


    @ODEmodel(
        x1'(t) = 3*x1(t) + 7*x2(t) - 6,
        x2'(t) = 9*x1(t) + 17*x2(t) - 16,
        y(t) = x2(t) + x1(t) + 1
    ), # [t1, -t1 - 1]

    @ODEmodel(
        x1'(t) = 3*x1(t) + 7*x2(t) - 6,
        x2'(t) = 9*x1(t) + 17*x2(t) - 16,
        y(t) = x2(t)^2 + x1(t) + 1
    ),

    @ODEmodel(
        x1'(t) = 2* x1(t) + x2(t) - 1,
        x2'(t) = x1(t) + x2(t) + 1,
        y(t) = x2(t) + x1(t) + 1
    ), # [t1, -t1 - 1]

    @ODEmodel(
        x1'(t) = 3*x1(t) + 7*x2(t) - 6,
        x2'(t) = 9*x1(t) + 17*x2(t) - 16,
        y(t) = x1(t)^2 + x1(t)^3 - x2(t)^2
    ), # [t1^2 - 1, t1 * (t1^2 - 1)]



    @ODEmodel(
        x1'(t) = 3*x1(t) + 7*x2(t) - 6,
        x2'(t) = 9*x1(t) + 17*x2(t) - 16,
        y(t) = x1(t)^2 + x2(t)^2 - 1
    ), #[ (1 - t1^2) // (1 + t1^2), (2 * t1) // ( 1 + t1^2)]

    @ODEmodel(
        x1'(t) = 3*x1(t) + 7*x2(t) - 6,
        x2'(t) = 9*x1(t) + 17*x2(t) - 16,
        y(t) = x1(t)*x2(t) + x2(t)^2 - x2(t) + 8 * x1(t) + 17
    ),

    @ODEmodel(
        x1'(t) = 3*x3(t) - x2(t) + x1(t) + 35 *  x3(t) ,
        x2'(t) = -3 * x2(t) + 2*(x1(t) - 35 * x2(t) - 6),
        x3'(t) = 6 * (x1(t) - 35 * x2(t) - 6)  +  32*x3(t),
        y(t) = x1(t)*x3(t) + x1(t)*x2(t) + x1(t) + 7 
    ),

    @ODEmodel(
        x1'(t) = 3*x1(t) - x2(t),
        x2'(t) = -3 * x2(t) + 2*x1(t) * x2(t),
        y(t) = x1(t)
    )

)

@testset "Testing against the standard algorithms" begin
    for c in cases
        @info c
        @test check_ansatz_modp(c, 2^31 - 1)
    end
end
