################################################################################
#
# benchmarks.jl
# =============
# 
# A benchmark file used to measure the timings of the computation of
# trisymmetric multiplication formulas in finite fields. 
#
################################################################################

using Nemo:FiniteField
using TriSym, Primes

path_trisym = dirname(dirname(pathof(TriSym)))
path_benchmark = joinpath(path_trisym, "benchmarks")
