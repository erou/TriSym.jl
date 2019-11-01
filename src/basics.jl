#######################################################
#
# Basic functions to handle/generate data
#
#######################################################

export multiplication_bilinear_map, make_conversion_dict

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
   
#######################################################
#
# Dictionnary
#
#######################################################

function make_conversion_dict(k::Nemo.FqNmodFiniteField)
    d = Dict{Nemo.fq_nmod, Nemo.nmod_mat}()
    for x in k
        d[deepcopy(x)] = trace_to_bil(x)
    end
    return d
end

#######################################################
#
# Multiplication bilinear map
#
#######################################################

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

function submul(a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)
    res = b*c
    ccall((:nmod_mat_sub, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}), res, a, res)
    return res
end

function submul!(res::Nemo.nmod_mat, a::Nemo.nmod_mat, b::Nemo.nmod_mat, c::UInt)
    ccall((:nmod_mat_scalar_mul, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, UInt), res, b, c)
    ccall((:nmod_mat_sub, :libflint), Nothing,
          (Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}, Ref{Nemo.nmod_mat}), res, a, res)
    return res
end
