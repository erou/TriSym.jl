module TriSym

using Nemo

import Base: +, *, copy, show, getindex, setindex!, parent, zero
import Nemo: ==, characteristic, base_ring, nmod, nmod_mat, fq_nmod, rank,
degree

#######################################################
#
# Types to represent data
#
#######################################################

include("types.jl")

#######################################################
#
# Basic functions to handle/generate data
#
#######################################################

include("basics.jl")

#######################################################
#
# Functions searching for decompositions
#
#######################################################

include("search.jl")
#include("search2.jl")

end
