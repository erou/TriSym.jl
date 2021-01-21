#######################################################
#
# Basic functions to handle/generate data
#
#######################################################

export multiplication_bilinear_map, make_conversion_dict, make_conversion_dict2

"""
    trace_to_bil(x::Nemo.fq_nmod)

Compute the bilinear map associated with the trace of an element.
"""
function trace_to_bil(x::Nemo.fq_nmod)
    k = parent(x)
    g = gen(k)
    n = degree(k)
    p::Int = characteristic(k)
    S = MatrixSpace(ResidueRing(ZZ, p), n, 1)
    a = S()
    for j in 1:n
        a[j, 1] = coeff(tr(x * g^(j-1)), 0)
    end
    return a * transpose(a)
end

function vector(x::Nemo.fq_nmod, S::Nemo.NmodMatSpace)
    a = S()
    for j in 1:nrows(S)
        a[j, 1] = coeff(x, j-1)
    end
    return a
end

function make_traces_matrix(k::Nemo.FqNmodFiniteField)
    n = degree(k)
    p::Int = characteristic(k)
    S = MatrixSpace(ResidueRing(ZZ, p), n, n)
    M = S()
    g = gen(k)
    for j in 1:n
        for i in 1:n
            M[i, j] = coeff(tr(g^(i+j-2)), 0)
        end
    end
    return M
end

function trace_to_bil(x::Nemo.fq_nmod, M::Nemo.nmod_mat, S::Nemo.NmodMatSpace)
    v = vector(x, S)
    a = M*v
    return a * transpose(a)
end

#######################################################
#
# Coeff manipulation
#
#######################################################

function coeff!(x::fq_nmod, j::Int, c::UInt)
    ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
          (Ref{fq_nmod}, Int, UInt), x, j, c)
end

coeff!(x::fq_nmod, j::Int, c::Int) = coeff!(x, j, UInt(c))

#######################################################
#
# Iterators over finite field
#
#######################################################

function Base.iterate(k::Nemo.FqNmodFiniteField)
    z = k()
    return (z, z)
end

function Base.iterate(k::Nemo.FqNmodFiniteField, state::Nemo.fq_nmod)
    p::UInt = characteristic(k)
    n = degree(k)
    l = n-1
    c0 = coeff(state, 0)
    if c0 < p-1
        coeff!(state, 0, c0+1)
    else
        for i in 0:n-1
            if coeff(state, i) != p-1
                l = i
                break
            end
        end

        for i in 0:l-1
            coeff!(state, i, 0)
        end
        coeff!(state, l, coeff(state, l)+1)
    end
    return state == 0 ? nothing : (state, state)
end

# iterate through elements without actually constructing the array
# Seems like a **very** poor idea though

struct AffineFieldElements
    parent::Nemo.FqNmodFiniteField
    coord::Int
    first::Nemo.fq_nmod
end

function AffineFieldElements(k::Nemo.FqNmodFiniteField, c::Int)
    z = k()
    coeff!(z, c-1, 1)
    return AffineFieldElements(k, c, z)
end

function Base.iterate(k::AffineFieldElements)
    z = k.first
    if z != 0
        zz = deepcopy(z)
        return (zz, zz)
    end
end

function Base.iterate(k::AffineFieldElements, state::Nemo.fq_nmod)
    p::UInt = characteristic(k.parent)
    n = degree(k.parent)
    l = n-1
    coord = k.coord
    g = k.parent()
    coeff!(g, coord-1, 1)
    if coord < n
        c0 = coeff(state, coord)
        if c0 < p-1
            coeff!(state, coord, c0+1)
        else
            for i in coord:n-1
                if coeff(state, i) != p-1
                    l = i
                    break
                end
            end

            for i in coord:l-1
                coeff!(state, i, 0)
            end
            coeff!(state, l, coeff(state, l)+1)
        end
        return state == g ? nothing : (state, state)
    end
end

function next(k::AffineFieldElements, state::Nemo.fq_nmod)
    nex = Base.iterate(k, deepcopy(state))
    if nex == nothing
        return k.parent()
    else
        return nex[1]
    end
end

function AffineFieldElements2(k::Nemo.FqNmodFiniteField, c::Int)
    p::BigInt = characteristic(k)
    d = degree(k)
    A = Array{Nemo.fq_nmod, 1}(undef, p^(d-c))
    j = 1
    for x in AffineFieldElements(k, c)
        A[j] = deepcopy(x)
        j += 1
    end
    return A
end

function Base.iterate(A::AbstractCommutativeAlgebra{N}) where N
    z = A()
    return (z, z)
end

function Base.iterate(A::AbstractCommutativeAlgebra{N},
                      state::AbstractCommutativeAlgebraElement{N}) where N
    p = characteristic(A)
    l = N
    c1 = state[1]
    if data(c1) < p-1
        state[1] += 1
    else
        for i in 1:N
            if state[i] != p-1
                l = i
                break
            end
        end

        for i in 1:l-1
            state[i] = 0
        end
        state[l] += 1
    end
    return state == 0 ? nothing : (state, state)
end

function AffineFieldElements2(A::AbstractCommutativeAlgebra{N}, c::Int) where N
    p::BigInt = characteristic(A)
    L = Array{nmod_mat, 1}(undef, p^(N-c))
    j = 1
    for x in A
        cont = false
        for i in 1:c-1
            if x[i] != 0
                cont = true
                break
            end
        end
        if cont
            continue
        elseif x[c] == 1
            L[j] = deepcopy(vector(x))
            j += 1
        end
    end
    return L
end
   
#######################################################
#
# Dictionnary
#
#######################################################

"""
    make_conversion_dict(k::Nemo.FqNmodFiniteField)

Compute a dictionnary mapping the elements of a finite field to their associated
bilinear map.
"""
function make_conversion_dict(k::Nemo.FqNmodFiniteField)
    d = Dict{Nemo.fq_nmod, Nemo.nmod_mat}()
    for x in k
        d[deepcopy(x)] = trace_to_bil(x)
    end
    return d
end

function make_conversion_dict2(k::Nemo.FqNmodFiniteField)
    d = Dict{Nemo.fq_nmod, Nemo.nmod_mat}()
    n = degree(k)
    p::Int = characteristic(k)
    T = make_traces_matrix(k)
    S = MatrixSpace(ResidueRing(ZZ, p), n, 1)
    for x in k
        d[deepcopy(x)] = trace_to_bil(x, T, S)
    end
    return d
end

"""
    make_conversion_dict(A::AbstractCommutativeAlgebra{N}, T::nmod_mat) where N

Compute a dictionnary mapping the elements of an algebra `A` to their associated
map.
"""
function make_conversion_dict(A::AbstractCommutativeAlgebra{N}, T::nmod_mat) where N
    M = multiplication_tensor(A)
    L =  T*M
    S = MatrixSpace(base_ring(M), N, N)()
    for i in 1:N
        for j in 1:i-1
            S[i, j] = S[j, i]
        end
        for j in i:N
            S[i, j] = L[1, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1]
        end
    end

    d = Dict{nmod_mat, nmod_mat}()

    for x in A
        v = deepcopy(vector(x))
        Sv = S*v
        d[v] = Sv*transpose(Sv)
    end

    return d
end

#######################################################
#
# Multiplication bilinear map
#
#######################################################

"""
    multiplication_bilinear_map(k::Nemo.FqNmodFiniteField)

Compute the bilinear map associated with the multiplication 
of a finite field `k`.
"""
function multiplication_bilinear_map(k::Nemo.FqNmodFiniteField)
    P = modulus(k) 
    n = degree(k)
    p::Int = characteristic(k)
    S = MatrixSpace(ResidueRing(ZZ, p), n, n)
    T = Array{Nemo.nmod_mat, 1}(undef, 2*n-1)
    for j in 1:2*n-1
        T[j] = S()
    end
    # This is **NOT** the same as `T = fill(S(), 2*n-1)`

    for j in 0:n-1
        for i in 0:n-1
            T[1+i+j][j+1, i+1] += 1
        end
    end

    l = 2*n-1
    while l > n
        for j in 0:n
            T[l-n+j] -= coeff(P, j)*T[l]
        end
        l -= 1
        pop!(T)
    end
    return BilinearMap(T, k)
end

"""
    multiplication_bilinear_map(A::AbstractCommutativeAlgebra{N}) where N

Compute the bilinear map assciated with the multiplication in the algebra `A`.
"""
function multiplication_bilinear_map(A::AbstractCommutativeAlgebra{N}) where N
    M = multiplication_tensor(A)
    S = MatrixSpace(base_ring(M), N, N)
    T = Array{nmod_mat, 1}(undef, N)
    for j in 1:N
        T[j] = S()
    end
    # This is **NOT** the same as `T = fill(S(), 2*n-1)`

    for i in 1:N
        XY = S()
        XY[i, i] = 1
        EE = M[:, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + 1]
        for k in 1:N
            T[k] += EE[k, 1] * XY
        end

        for j in i+1:N
            XY = S()
            XY[i, j], XY[j, i] = 1, 1
            EE = M[:, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1]
            for k in 1:N
                T[k] += EE[k, 1] * XY
            end
        end
    end

    return BilinearMap2(T, A)
end

#######################################################
#
# Ad hoc function to save allocations
#
#######################################################

"""
    submul(a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)

Compute `a - b*c` with 1 allocation instead of 2.
"""
function submul(a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)
    res = b*c
    ccall((:nmod_mat_sub, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}), res, a, res)
    return res
end

"""
    submul!(res::Nemo.nmod_mat, a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)

Compute `a - b*c` and put the result in `res`.

Note
====

No allocations performed.
"""
function submul!(res::Nemo.nmod_mat, a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)
    ccall((:nmod_mat_scalar_mul, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, UInt), res, b, c)
    ccall((:nmod_mat_sub, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}), res, a, res)
end

"""
    addmul(a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)

Compute `a + b*c` with 1 allocation instead of 2.
"""
function addmul(a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)
    res = b*c
    ccall((:nmod_mat_add, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}), res, a, res)
    return res
end

"""
    addmul!(res::Nemo.nmod_mat, a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)

Compute `a + b*c` and put the result in `res`.

Note
====

No allocations performed.
"""
function addmul!(res::Nemo.nmod_mat, a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)
    ccall((:nmod_mat_scalar_mul, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, UInt), res, b, c)
    ccall((:nmod_mat_add, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}), res, a, res)
end

#######################################################
#
# Dictionnary for time / memory trade off
#
#######################################################

"""
    _make_small_trisym_coord(elems::Array{Nemo.fq_nmod, 1}, count::Int,
    res::Dict{Nemo.nmod_mat, Array{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}, 1}},
    current_sum::Nemo.nmod_mat, decomp::Array{Tuple{Nemo.fq_nmod, Int}, 1},
    dic::Dict{Nemo.fq_nmod, Nemo.nmod_mat})

Compute a dictionnary maping bilinear maps to their trysymmetric decompositions.

Note
====

This is a private function that should not be called. See the function
`make_small_trisym_coord1`.
"""
function _make_small_trisym_coord(elems::Array{Nemo.fq_nmod, 1}, count::Int,
         res::Dict{Nemo.nmod_mat, Array{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}, 1}},
         current_sum::Nemo.nmod_mat, decomp::Array{Tuple{Nemo.fq_nmod, Int}, 1},
         dic::Dict{Nemo.fq_nmod, Nemo.nmod_mat})

    p::UInt = characteristic(parent(first(keys(dic))))
    if count == 1
        for y in elems
            for c in 1:p-1
                new_decomp = Base.copy(decomp)
                push!(new_decomp, (y, c))
                #new_sum = current_sum + c*dic[y] # 2 allocations here, only one needed
                new_sum = addmul(current_sum, dic[y], c) # 1 allocation here
                if haskey(res, new_sum)
                    push!(res[new_sum], Tuple(new_decomp))
                else
                    res[new_sum] = Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}[Tuple(new_decomp)]
                end
            end
        end
    else#  count > 1
        for i in 1:length(elems)
            y = elems[i]
            for c in 1:p-1
                new_decomp = Base.copy(decomp)
                push!(new_decomp, (y, c))
                #new_sum = current_sum + c*dic[y] # 2 allocations here, only one needed
                new_sum = addmul(current_sum, dic[y], c) # 1 allocation here
                if haskey(res, new_sum)
                    push!(res[new_sum], Tuple(new_decomp))
                else
                    res[new_sum] = Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}[Tuple(new_decomp)]
                end
                _make_small_trisym_coord(elems[i+1:end], count-1, res,
                                         new_sum, new_decomp, dic)
            end
        end
    end
end

"""
    make_small_trisym_coord(k::Nemo.FqNmodFiniteField, rank::Int,
                            dic::Dict{Nemo.fq_nmod, Nemo.nmod_mat}, j::Int)

Compute a dictionnary maping bilinear maps to their trysymmetric decompositions.

Parameters
==========

- `rank` the length of the computed decompositions
- `k` the finite fields considered
- `j` the coordinate
- `dic` the dictionnary mapping the elements of `k` to their associated bilinear
  maps

Note
====

Only small trisymmetric decompositions are computed, because of the exponential
complexity of doing so.
"""
function make_small_trisym_coord(k::Nemo.FqNmodFiniteField, rank::Int,
                                 dic::Dict{Nemo.fq_nmod, Nemo.nmod_mat}, j::Int)

    elems = AffineFieldElements2(k, j)
    S = parent(dic[gen(k)])
    current_sum = S()
    res_tmp = Dict{Nemo.nmod_mat, Array{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}, 1}}()
    decomp = Tuple{Nemo.fq_nmod, Int}[]

    _make_small_trisym_coord(elems, rank, res_tmp, current_sum, decomp, dic)

    res = Dict{Nemo.nmod_mat, Tuple{Vararg{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}}}}()

    for x in keys(res_tmp)
        res[x] = Tuple(res_tmp[x])
    end

    res[S()] = ((),)

    return res
end

#######################################################
#
# sort through fq_nmod arrays
#
#######################################################

function Base.isless(a::Nemo.fq_nmod, b::Nemo.fq_nmod)
    @assert parent(a) == parent(b)
    n = degree(parent(a))
    for j in 0:n-1
        if coeff(a, j) < coeff(b, j)
            return true
        elseif coeff(a, j) != coeff(b, j)
            return false
        end
    end
    return false
end

#######################################################
#
# Basic functions to handle AbstractCommutativeAlgebra
#
#######################################################

function (A::AbstractCommutativeAlgebra{N})() where N
    return AbstractCommutativeAlgebraElement(A)
end

function (A::AbstractCommutativeAlgebra{N})(L::NTuple{N, Int}) where N
    a = A()
    for i in 1:N
        a[i] = L[i]
    end
    return a
end

function +(a::AbstractCommutativeAlgebraElement{N},
            b::AbstractCommutativeAlgebraElement{N}) where N

    @assert parent(a) == parent(b)
    A = parent(a)

    return AbstractCommutativeAlgebraElement(vector(a)+vector(b), A)
end

function +(a::AbstractCommutativeAlgebraElement{N}, b::nmod_mat) where N

    A = parent(a)
    return AbstractCommutativeAlgebraElement(vector(a)+b, A)
end

function *(a::AbstractCommutativeAlgebraElement{N},
           b::AbstractCommutativeAlgebraElement{N}) where N

    @assert parent(a) == parent(b)
    A = parent(a)
    M = multiplication_tensor(A)
    c = A()
    for i in 1:N
        c += a[i]*b[i]*M[:, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + 1]
        for j in i+1:N
            c += (a[i]*b[j]+a[j]*b[i])*M[:, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1]
        end
    end
    return c
end

#= ######################################################
#
# Test functions
#   or
#     gotta go fast!
#
#######################################################

function test1(A::Array{nmod_mat, 1}, d::Dict{nmod_mat, nmod_mat})
    for x in A
        d[x]
    end
end

function test2(A::Array{NTuple{4, UInt8}, 1}, d::Dict{NTuple{4, UInt8}, nmod_mat})
    for x in A
        d[x]
    end
end

function test3(A::Array{NTuple{4, UInt8}, 1}, d::NTuple{3^4, nmod_mat})
    for x in A
        d[x]
    end
end

function test4(A::Array{UInt32, 1}, d::NTuple{3^4, nmod_mat})
    for x in A
        d[x]
    end
end
=#

#######################################################
#
# Basic functions to handle CommutativeAlgebra
#
# those were not so fast in the end...
#
#######################################################

function Base.iterate(A::CommutativeAlgebra{N, U}) where {N, U <: Integer}

    z = MatrixSpace(ResidueRing(ZZ, Int(characteristic(A))), N, 1)()
    return (z, z)
end

function Base.iterate(A::CommutativeAlgebra{N}, state::nmod_mat) where {N, U <: Integer}
    p = characteristic(A)
    l::Int = N
    c1 = state[1, 1]
    if data(c1) < p-1
        state[1, 1] += 1
    else
        for i in 1:N
            if state[i, 1] != p-1
                l = i
                break
            end
        end

        for i in 1:l-1
            state[i, 1] = 0
        end
        state[l, 1] += 1
    end
    return state == 0 ? nothing : (state, state)
end

function make_conversion_tuple(A::CommutativeAlgebra{N, U}, T::nmod_mat) where {N, U <: Integer}
    M = multiplication_tensor(A)
    p = characteristic(A)
    L =  T*M
    S = MatrixSpace(base_ring(M), N, N)()
    for i in 1:N
        for j in 1:i-1
            S[i, j] = S[j, i]
        end
        for j in i:N
            S[i, j] = L[1, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1]
        end
    end

    d = Array{nmod_mat, 1}(undef, p^N)

    k = 1
    for x in A
        v = deepcopy(x)
        Sv = S*v
        d[k] = Sv*transpose(Sv)
        k += 1
    end

    return Tuple(d)
end

function AffineElements(A::CommutativeAlgebra{N, U}, c::U) where {N, U <: Integer}
    p = characteristic(A)
    o = one(U)
    t = p^(N-c)
    F = Array(zero(U):t-o)
    Ans = @. p^(c-o) + p^c * F
    return Tuple(Ans)
end

function multiplication_bilinear_map(A::CommutativeAlgebra{N, U}) where {N, U <: Integer}
    M = multiplication_tensor(A)
    S = MatrixSpace(base_ring(M), N, N)
    T = Array{nmod_mat, 1}(undef, N)
    for j in 1:N
        T[j] = S()
    end
    # This is **NOT** the same as `T = fill(S(), 2*n-1)`

    for i in 1:N
        XY = S()
        XY[i, i] = 1
        EE = M[:, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + 1]
        for k in 1:N
            T[k] += EE[k, 1] * XY
        end

        for j in i+1:N
            XY = S()
            XY[i, j], XY[j, i] = 1, 1
            EE = M[:, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1]
            for k in 1:N
                T[k] += EE[k, 1] * XY
            end
        end
    end

    return BilinearMap3(T, A)
end

function translate(x::U, j::Int, N::Int, p::U) where U <: Integer
    y = div(x, p^j)
    t = Array{U, 1}(undef, N-j)

    for i in 1:N-j
        t[i] = y % p
        y ÷= p 
    end
    return Tuple(t)
end

#######################################################
#
# sort through nmod_mat arrays
#
#######################################################

function Base.isless(a::nmod_mat, b::nmod_mat)
    @assert parent(a) == parent(b)
    n, m = nrows(a), ncols(a)
    for j in 1:m
        for i in 1:n
            if data(a[i, j]) < data(b[i, j])
                return true
            elseif data(a[i, j]) != data(b[i, j])
                return false
            end
        end
    end
    return false
end

#######################################################
#
# Functions for a vector space strategy
#
#######################################################

"""
    make_conversion_dict_vectors(A::AbstractCommutativeAlgebra{N}, T::nmod_mat) where N

Compute a dictionnary mapping the elements of the algebra `A` with their
associated bilinear map, via the residue `T`.
"""
function make_conversion_dict_vectors(A::AbstractCommutativeAlgebra{N}, T::nmod_mat) where N
    M = multiplication_tensor(A)
    L =  T*M
    S = MatrixSpace(base_ring(M), N, N)()
    for i in 1:N
        for j in 1:i-1
            S[i, j] = S[j, i]
        end
        for j in i:N
            S[i, j] = L[1, div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1]
        end
    end

    d = Dict{nmod_mat, nmod_mat}()

    for x in A
        v = deepcopy(vector(x))
        Sv = S*v
        d[v] = kronecker_product(v, kronecker_product(Sv, Sv))
    end

    return d
end

"""
    multi_rank(x::nmod_mat, N::Int)

Compute the ranks of the `N` matrices of size N×N encoded in the vector `x`
that is of size `N^3`.
"""
function multi_rank(x::nmod_mat, N::Int)

    S = MatrixSpace(base_ring(parent(x)), N, N)
    ranks = Array{Int, 1}(undef, N)
    for j in 1:N
        A = S()
        for i in 1:N
            for k in 1:N
                A[k, i] = x[(j-1)*N^2 + (i-1)*N + k, 1]
            end
        end
        ranks[j] = rank(A)
    end
    return ranks
end

"""
    rank(X::Array{nmod_mat, 1})

Compute the rank of the matrix whose columns are the vectors in `X`.
"""
function rank(X::Array{nmod_mat, 1})
    l = length(X)
    if l > 0
        n = nrows(X[1])
        S = MatrixSpace(base_ring(X[1]), n, l)()
        for j in 1:l
            x = X[j]
            for i in 1:n
                S[i, j] = x[i, 1]
            end
        end

        return rank(S)

    else
        return 0
    end
end

"""
    normalize!(x::nmod_mat)

Find the first nonzero coordinate of `x` and divide `x` by the value of this
coordinate. Return the index of this coordinate.

Note
====

The element `x` is supposed to be a column matrix.
"""
function normalize!(x::nmod_mat)

    # Find the first nonzero entry
    j = 0
    n = nrows(x)
    for k in 1:n
        if x[k, 1] != 0
            j = k
            break
        end
    end

    # Divide the other entries by its value
    if j != 0 
        x0 = x[j, 1]^(-1)
        for k in j:n
            x[k, 1] *= x0
        end
    end

    return j
end

"""
    reduce!(x::nmod_mat, W::nmod_mat)

Reduce the vector `x` modulo the vector space represented by `W`.

Note
====

W must be in **reduced column echelon form** and *square**.
"""
function reduce!(x::nmod_mat, W::nmod_mat)
    n = nrows(x)

    for j in 1:n
        if x[j, 1] != 0 && W[j, j] != 0
            y = x[j, 1]
            for i in j:n
                x[i, 1] -= y*W[i, j]
            end
        end
    end
end

"""
    reduce_column!(x::nmod_mat, y::nmod_mat)
Reduce the vector column matrix `x` by the column matrix `y`.

Note
====

The matrix `y` is assumed to be normalized.
"""
function reduce_column!(x::nmod_mat, y::nmod_mat)

    # Find the first nonzero entry of `y`
    j = 0
    n = nrows(y)
    for k in 1:n
        if y[k, 1] != 0
            j = k
            break
        end
    end

    # If the corresponding entry in `x` is nonzero, substract in order to
    # nullify it
    x0 = x[j, 1]
    if x0 != 0
        for k in j:n
            x[k, 1] -= x0 * y[k, 1]
        end
    end
end

"""
    reduce!(X::Array{nmod_mat, 1}, W::nmod_mat)

Reduce each element in the array `X` modulo the vector space represented
by `W`. Then delete the doubles.

Note
====

W must be in **reduced column echelon form**.
"""
function reduce!(X::Array{nmod_mat, 1}, W::nmod_mat)

    # Reduction of all entries
    for x in X
        reduce!(x, W)
    end

    # We sort the result
    sort!(X)

    # Then we delete the doubles
    j = 1
    l = length(X)
    while j < l
        if X[j] == X[j+1]
            deleteat!(X, j+1)
            l -= 1
        else
            j += 1
        end
    end

    # Then we remove the zeros, if there are any
    if l > 0
        if X[1] == 0
            deleteat!(X, 1)
        end
    end
end

"""
    reduce_column!(X::Array{nmod_mat, 1}, y::nmod_mat)

Reduce each element in the array `X` modulo the vector space represented
by `W`. Then delete the doubles.

Note
====

W must be in **reduced column echelon form**.
"""
function reduce_column!(X::Array{nmod_mat, 1}, y::nmod_mat)

    # Reduction of all entries
    for x in X
        reduce_column!(x, y)
    end

    # We sort the result
    sort!(X)

    # Then we delete the doubles
    j = 1
    l = length(X)
    while j < l
        if X[j] == X[j+1]
            deleteat!(X, j+1)
            l -= 1
        else
            j += 1
        end
    end

    # Then we remove the zeros, if there are any
    if l > 0
        if X[1] == 0
            deleteat!(X, 1)
        end
    end
end

"""
    add_vector!(W::nmod_mat, x::nmod_mat)

Add the vector `x` in `W`, respecting the reduced column echelon form.
"""
function add_vector!(W::nmod_mat, x::nmod_mat)

    # We reduce x by W
    reduce!(x, W)

    # We normalize x
    j = normalize!(x)

    # If j != 0 it means that x is not in W thus we should add it
    if j != 0
        n = nrows(x)
        for i in j:n
            W[i, j] = x[i, 1]
        end

        # We should also reduce W by x now
        for i in 1:j-1
            if W[j, i] != 0
                w0 = W[j, i]
                for k in j:n
                    W[k, i] -= w0 * W[k, j]
                end
            end
        end
    end
end

"""
    construct_matrix(d::Dict{nmod_mat, nmod_mat},
                     A::AbstractCommutativeAlgebra{N}) where N

Construct the matrix which columns are the vectors representing the normalized
trisymmetric rank one tensors.
"""
function construct_matrix(d::Dict{nmod_mat, nmod_mat},
                          A::AbstractCommutativeAlgebra{N}) where N

    Z = base_ring(parent(first(keys(d))))
    p::Int = Z.n
    S = MatrixSpace(Z, N^3, div(p^N-1, p-1))()
    cpt = 1
    for j in 1:N
        G = AffineFieldElements2(A, j)
        for x in G
            v = d[x]
            for k in 1:N^3
                S[k, cpt] = v[k, 1]
            end
            cpt += 1
        end
    end
    return S
end

"""
    intersection(W::nmod_mat, G::nmod_mat)

Compute an array which components are vectors that are in the intersection of
`G` and the vector space represented by `W`.
"""
function intersection(W::nmod_mat, G::nmod_mat)
    I = W*G
    L = nmod_mat[]
    for j in 1:ncols(G)
        colj = G[:, j]
        if I[:, j] == colj
            push!(L, colj)
        end
    end
    return L
end

"""
    multiplication_vector(A::AbstractCommutativeAlgebra{N}) where N

Compute the vector representing the tensor of multiplication of the algebra `A`.
"""
function multiplication_vector(A::AbstractCommutativeAlgebra{N}) where N
    B = multiplication_bilinear_map(A)
    Z = base_ring(multiplication_tensor(A))
    S = MatrixSpace(Z, N^3, N^3)()
    cpt = 1
    for j in 1:N
        for i in 1:N
            for k in 1:N
                S[cpt, 1] = B[j][k, i]
                cpt += 1
            end
        end
    end
    return S
end

"""
    mycopy(X::Array{nmod_mat, 1})

Return an array of elements which are copies of the elements in `X`.

Note
====

This is the same as `copy`, and probably not the same as `deepcopy`.
"""
function mycopy(X::Array{nmod_mat, 1})
    Y = similar(X)
    for (j, x) in enumerate(X)
        Y[j] = copy(x)
    end
    return Y
end

"""
    rank_one_basis(W::nmod_mat, G::nmod_mat)

Compute a basis of `W` composed of the columns of `G`.
"""
function rank_one_basis(W::nmod_mat, G::nmod_mat)
    I = intersection(W, G)
    r = rank(W)
    S = similar(W)
    add_vector!(S, copy(I[1]))
    base = Array{nmod_mat, 1}(undef, r)
    base[1] = I[1]
    j = 2
    k = 2
    while j <= r
        y = copy(I[k])
        reduce!(y, S)
        if y != 0
            add_vector!(S, y)
            base[j] = I[k]
            j += 1
        end
        k += 1
    end
    return base
end

"""
    coordinates(x::nmod_mat, basis::Array{nmod_mat, 1})

Compute the coordinates of the vector `x` in the basis `basis`.
"""
function coordinates(x::nmod_mat, basis::Array{nmod_mat, 1})
    n = nrows(x)
    m = length(basis)+1
    S = MatrixSpace(base_ring(x), n, m)()
    for j in 1:m-1
        b = basis[j]
        for i in 1:n
            S[i, j] = b[i, 1]
        end
    end

    for i in 1:n
        S[i, m] = x[i, 1]
    end

    rref!(S)

    return S[:, m]
end

function get_formula(w::nmod_mat, W::nmod_mat, G::nmod_mat, d::Dict{nmod_mat,
                     nmod_mat})
    basis = rank_one_basis(W, G)
    coord = coordinates(w, basis)

    ke = [v for v in keys(d)]
    va = [w for w in values(d)]
    
    return [(ke[indexin([basis[j]], va)[1]], coord[j, 1]) for j in 1:length(basis)] 
end

#######################################################
#
# Type tryouts
#
#######################################################

function Base.iterate(A::ComAlgebra)
    z = A()
    return (z, z)
end

function Base.iterate(A::ComAlgebra, state::nmod_mat)
    p = characteristic(A)
    N = degree(A)
    l = N
    c1 = state[1, 1]
    if data(c1) < p-1
        state[1, 1] += 1
    else
        for i in 1:N
            if state[i, 1] != p-1
                l = i
                break
            end
        end

        for i in 1:l-1
            state[i, 1] = 0
        end
        state[l, 1] += 1
    end
    return state == 0 ? nothing : (state, state)
end

function AffineFieldElements2(A::ComAlgebra, c::Int)
    p::BigInt = characteristic(A)
    N = degree(A)
    L = Array{nmod_mat, 1}(undef, p^(N-c))
    j = 1
    for x in A
        cont = false
        for i in 1:c-1
            if x[i, 1] != 0
                cont = true
                break
            end
        end
        if cont
            continue
        elseif x[c, 1] == 1
            L[j] = deepcopy(x)
            j += 1
        end
    end
    return L
end
