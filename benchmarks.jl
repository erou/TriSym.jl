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

function small_margin(N, n, bound, margin, str = "fixed-margin")

    path = joinpath(path_benchmark, str)
    for i in 1:n
        path = joinpath(path, string("-", margin[i]))
    end
    path = joinpath(path, ".txt")
    io = open(path, "w+")

    P = primes(3, N)
    for p in P
        k, x = FiniteField(p, n, "x")
        B = multiplication_bilinear_map(k)
        res = @timed make_conversion_dict(k)
        d, t1 = res[1], res[2]
        L = Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}[]
        t2 = @elapsed tri_symmetric_search(B, d, L, margin, bound)
        write(io, p, ",", t1, ",", t2, "\n")
    end

    close(io)

end
