include("solver.jl")

using Test
using StructuralIdentifiability: _reduce_mod_p, str_to_var, eval_at_dict

@testset "Testing against the standard algorithms" begin
    ode = @ODEmodel(
        x1'(t) = x1(t)^2 - 3 * x1(t) * x2(t) + 2 * x2(t)^2 - 3 * x1(t) + 2 * x2(t) - 7,
        x2'(t) = x1(t)^2 + x1(t) * x2(t) - x2(t)^2 + x1(t) + 1,
        y(t) = x1(t)
    )
    
    p = 2^30 + 3
    
    io_tocheck = solve_with_love_and_support(ode, p)
    
    # computing the IO-equation by the standard algorithm
    io_correct = first(values(find_ioequations(ode)))
    io_correct = _reduce_mod_p(io_correct, p)
    
    # bringing them to the same ring
    yvars = [str_to_var(s, parent(io_correct)) for s in ("y(t)_0", "y(t)_1", "y(t)_2")]
    xvars = [str_to_var(s, parent(io_tocheck)) for s in ("x1", "x1'", "x1''")]
    io_correct = eval_at_dict(io_correct, Dict(yvars .=> xvars))
    
    # checking that they are proportional
    quot, rem = divrem(io_tocheck, io_correct)
    @test iszero(rem)
    @test iszero(total_degree(quot))
end
