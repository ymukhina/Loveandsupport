module DiffMinPoly

include("reduce_and_reorder.jl")
include("solver_love_and_support.jl")
export eliminate, eliminate_with_love_and_support, eliminate_with_love_and_support_modp, rand_ode

end # module DiffMinPoly
