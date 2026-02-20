module DiffMinPoly

include("epsilon.jl")
include("matrix.jl")
include("utils.jl")
include("helpers.jl")
include("solver_love_and_support.jl")
include("../test/test_utils.jl")


export eliminate, eliminate_with_love_and_support, eliminate_with_love_and_support_modp, rand_ode, Epsilon, rand_ode_x

end # module DiffMinPoly
