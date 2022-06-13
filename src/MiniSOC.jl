module MiniSOC

using LinearAlgebra

export heuristic, GreedyPowertwo, check_SOC, traversal, bruteforce, config

include("heuristic.jl")
include("traversal.jl")
include("bruteforce.jl")

end
