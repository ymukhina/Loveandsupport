module DiffMinPoly

include("matrix_init.jl")
include("utils.jl")
include("epsilon.jl")
include("solver_love_and_support.jl")
include("../test/test_utils.jl")

export eliminate, eliminate_with_love_and_support, eliminate_with_love_and_support_modp, rand_ode, Epsilon, rand_ode_x

end # module DiffMinPoly
