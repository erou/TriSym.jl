module TriSym

using Nemo

import Nemo: ==

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

end
