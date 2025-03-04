module DiffMinPoly

include("matrix_init.jl")
include("epsilon.jl")
include("solver_love_and_support.jl")
export eliminate, eliminate_with_love_and_support, eliminate_with_love_and_support_modp, rand_ode, Epsilon

end # module DiffMinPoly
