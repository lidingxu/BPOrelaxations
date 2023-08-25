module BPOrelaxations

const is_debug = false
const epsilon = 1e-6

using DataStructures
using JuMP
using MutableNamedTuples
using Graphs 
using GraphsFlows
using MosekTools

include("Utils.jl")
include("Reader.jl")
include("Problem.jl")
include("Separator.jl")
include("Submodular.jl")
include("Relaxations.jl")

export readmc, readBQP, relax, Option, parseOption, Stat, Problem
end # module BPOrelaxations
