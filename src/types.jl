##############################################################################
#
# Type to represent bilinear maps
# ===============================
#
# Presentation
# ============
#
# There are 3 different couples of types, each time coming with a type for the
# algebra and a type for the bilinear maps, that is often just a tuple with
# matrices and a reference to the algebra involved.
# 
# - BilinearMap : for bilinear maps between finite fields
# - BilinearMap2 : for bilinear maps between "abstract" algebras given by their
#   multiplication tensor
# - BilinearMap3 : for bilinear maps between finite fields, trying to use a
#   better representation that the first one
#
# Goal
# ====
#
# These functions provide the tools to manipulate bilinear maps, print them,
# etc.
#
##############################################################################

export BilinearMap, BilinearMap2, BilinearMap3

##############################################################################
#
# BilinearMap : the type for bilinear maps between Nemo's finite fields
# =====================================================================
#
##############################################################################

struct BilinearMap{N}
    coord::NTuple{N, Nemo.nmod_mat}
    base_ring::Nemo.FqNmodFiniteField

    function BilinearMap(k::Nemo.FqNmodFiniteField)
        N = degree(k)
        p::Int = characteristic(k)
        S = MatrixSpace(ResidueRing(ZZ, p), N, N)
        A = Array{Nemo.nmod_mat, 1}(undef, N)
        for j in 1:N
            A[j] = S()
        end
        return new{N}(Tuple(A), k)
    end

    function BilinearMap(A::Array{Nemo.nmod_mat, 1}, k::FqNmodFiniteField)
        return new{length(A)}(Tuple(A), k)
    end

    function BilinearMap(T::NTuple{N, Nemo.nmod_mat}, k::FqNmodFiniteField) where N
        return new{length(T)}(T, k)
    end
end

coord(b::BilinearMap{N}) where N = b.coord

base_ring(b::BilinearMap{N}) where N = b.base_ring

show(io::IO, b::BilinearMap{N}) where N = show(io, coord(b))

getindex(b::BilinearMap{N}, i::Int) where N = getindex(coord(b), i)

setindex!(b::BilinearMap{N}, i::Int) where N = getindex(coord(b), i)

zero(b::BilinearMap{N}) where N = BilinearMap(base_ring(b))

function check_rings(a::BilinearMap{N}, b::BilinearMap{N}) where N
    base_ring(a) != base_ring(b) && error("Bilinear maps do not have the same base ring")
end

function ==(a::BilinearMap{N}, b::BilinearMap{N}) where N
    check_rings(a, b)
    for j in 1:N
        if a[j] != b[j]
            return false
        end
    end
    return true
end

function iszero(b)
    a = zero(b)
    return a == b
end

function ==(a::BilinearMap{N}, i::Int) where N
    i != 0 && error("Cannot compare a bilinear map with a nonzero integer")
    return iszero(a)
end

function add!(b::BilinearMap{N}, decomp::Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}},
             j::Int, d::Dict{Nemo.fq_nmod, Nemo.nmod_mat}) where N

    zero!(b[j])

    for (x, c) in decomp
        M = -c*d[x]
        for i in j+1:N
            Nemo.add!(b[i], b[i], coeff(x, i-1)*M)
        end
    end
end

function copy(b::BilinearMap{N}) where N
    return BilinearMap(deepcopy(coord(b)), base_ring(b))
end

##############################################################################
#
# Type to represent abstract commutative algebras
# ===============================================
#
##############################################################################

struct AbstractCommutativeAlgebra{N}
    multiplication_tensor::Nemo.nmod_mat
    characteristic::Int

    function AbstractCommutativeAlgebra(p::Int, M::Nemo.nmod_mat)
        N = nrows(M)
        @assert ncols(M) == div(N*(N+1), 2)
        return new{N}(M, p)
    end

    function AbstractCommutativeAlgebra(k::FqNmodFiniteField)
        N = degree(k)
        p::Int = characteristic(k)
        M = MatrixSpace(ResidueRing(ZZ, p), N, div(N*(N+1), 2))()
        x = gen(k)
        for i in 1:N
            for j in i:N
                z = x^(i+j-2)
                col = div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1
                for l in 1:N
                    M[l, col] = coeff(z, l-1)
                end
            end
        end
        return new{N}(M, p)
    end
end

multiplication_tensor(A::AbstractCommutativeAlgebra{N}) where N = A.multiplication_tensor
characteristic(A::AbstractCommutativeAlgebra{N}) where N = A.characteristic

struct AbstractCommutativeAlgebraElement{N}
    vector::Nemo.nmod_mat
    parent::AbstractCommutativeAlgebra

    function AbstractCommutativeAlgebraElement(A::AbstractCommutativeAlgebra{N}) where N
        R = base_ring(multiplication_tensor(A))
        v = MatrixSpace(R, N, 1)()
        return new{N}(v, A)
    end

    function AbstractCommutativeAlgebraElement(v::Nemo.nmod_mat,
                                               A::AbstractCommutativeAlgebra{N}) where N
        return new{N}(v, A)
    end
end

vector(a::AbstractCommutativeAlgebraElement{N}) where N = a.vector
parent(a::AbstractCommutativeAlgebraElement{N}) where N = a.parent

show(io::IO, a::AbstractCommutativeAlgebraElement{N}) where N = show(io, vector(a))
getindex(a::AbstractCommutativeAlgebraElement, i::Int) where N = vector(a)[i, 1]
setindex!(a::AbstractCommutativeAlgebraElement, c::Int, i::Int) where N =
setindex!(vector(a), c, i, 1)
setindex!(a::AbstractCommutativeAlgebraElement, c::nmod, i::Int) where N =
setindex!(vector(a), c, i, 1)

zero(a::AbstractCommutativeAlgebraElement{N}) where N = parent(a)()

function ==(a::AbstractCommutativeAlgebraElement{N}, b::AbstractCommutativeAlgebraElement{N}) where N
    @assert parent(a) == parent(b)
    return vector(a) == vector(b)
end

function iszero(a::AbstractCommutativeAlgebraElement{N}) where N
    b = zero(a)
    return a == b
end

function ==(a::AbstractCommutativeAlgebraElement{N}, i::Int) where N
    i != 0 && error("Cannot compare an element with a nonzero integer")
    return iszero(a)
end

##############################################################################
#
# Type to represent bilinear maps between abstract algebras
# =========================================================
#
##############################################################################

struct BilinearMap2{N}
    coord::NTuple{N, nmod_mat}
    base_ring::AbstractCommutativeAlgebra{N}

    function BilinearMap2(A::AbstractCommutativeAlgebra{N}) where N
        S = MatrixSpace(base_ring(multiplication_tensor(A)), N, N)
        L = Array{nmod_mat, 1}(undef, N)
        for j in 1:N
            L[j] = S()
        end
        return new{N}(Tuple(L), A)
    end

    function BilinearMap2(L::Array{nmod_mat, 1},
                         A::AbstractCommutativeAlgebra{N}) where N
        return new{N}(Tuple(L), A)
    end

    function BilinearMap2(T::NTuple{N, nmod_mat}, A::AbstractCommutativeAlgebra{N}) where N
        return new{N}(T, A)
    end
end

coord(b::BilinearMap2{N}) where N = b.coord

base_ring(b::BilinearMap2{N}) where N = b.base_ring

show(io::IO, b::BilinearMap2{N}) where N = show(io, coord(b))

getindex(b::BilinearMap2{N}, i::Int) where N = getindex(coord(b), i)

zero(b::BilinearMap2{N}) where N = BilinearMap2(base_ring(b))

function check_rings(a::BilinearMap2{N}, b::BilinearMap2{N}) where N
    base_ring(a) != base_ring(b) && error("Bilinear maps do not have the same base ring")
end

function ==(a::BilinearMap2{N}, b::BilinearMap2{N}) where N
    check_rings(a, b)
    for j in 1:N
        if a[j] != b[j]
            return false
        end
    end
    return true
end

function ==(a::BilinearMap2{N}, i::Int) where N
    i != 0 && error("Cannot compare a bilinear map with a nonzero integer")
    return iszero(a)
end

function add!(b::BilinearMap2{N}, decomp::Tuple{Vararg{Tuple{nmod_mat, Int}}},
             j::Int, d::Dict{nmod_mat, Nemo.nmod_mat}) where N

    zero!(b[j])

    for (x, c) in decomp
        M = -c*d[x]
        for i in j+1:N
            Nemo.add!(b[i], b[i], x[i, 1]*M)
        end
    end
end

function copy(b::BilinearMap2{N}) where N
    return BilinearMap2(deepcopy(coord(b)), base_ring(b)) # copy sufficient?
end

##############################################################################
#
# Type to represent (optimised?) commutative algebras
# ===================================================
#
##############################################################################

function choose_type(p::Int, N::Int)
    l = log2(BigInt(p)^N)
    if l < 8
        return UInt8
    elseif l < 16
        return UInt16
    elseif l < 32
        return UInt32
    elseif l < 64
        return UInt64
    elseif l < 128
        return UInt128
    else
        error("Not implemented")
    end
end

struct CommutativeAlgebra{N, U <: Integer}
    multiplication_tensor::Nemo.nmod_mat
    characteristic::U

    function CommutativeAlgebra(p::Int, M::Nemo.nmod_mat)
        N = nrows(M)
        U = choose_type(p, N)
        return new{N, U}(M, U(p))
    end

    function CommutativeAlgebra(k::FqNmodFiniteField)
        N = degree(k)
        p::Int = characteristic(k)
        U = choose_type(p, N)
        M = MatrixSpace(ResidueRing(ZZ, p), N, div(N*(N+1), 2))()
        x = gen(k)
        for i in 1:N
            for j in i:N
                z = x^(i+j-2)
                col = div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1
                for l in 1:N
                    M[l, col] = coeff(z, l-1)
                end
            end
        end
        return new{N, U}(M, U(p))
    end
end

multiplication_tensor(A::CommutativeAlgebra{N, U}) where {N, U <: Integer} = A.multiplication_tensor
characteristic(A::CommutativeAlgebra{N, U}) where {N, U <: Integer} = A.characteristic

##############################################################################
#
# Type to represent bilinear maps between ad hoc finite fields
# ============================================================
#
##############################################################################

struct BilinearMap3{N, U <: Integer}
    coord::NTuple{N, nmod_mat}
    base_ring::CommutativeAlgebra{N, U}

    function BilinearMap3(A::CommutativeAlgebra{N, U}) where {N, U <: Integer}
        S = MatrixSpace(base_ring(multiplication_tensor(A)), Int(N), Int(N))
        L = Array{nmod_mat, 1}(undef, N)
        for j in 1:N
            L[j] = S()
        end
        return new{N, U}(Tuple(L), A)
    end

    function BilinearMap3(L::Array{nmod_mat, 1},
                          A::CommutativeAlgebra{N, U}) where {N, U <: Integer}
        return new{N, U}(Tuple(L), A)
    end

    function BilinearMap3(T::NTuple{N, nmod_mat}, A::CommutativeAlgebra{N, U}) where {N, U <: Integer}
        return new{N, U}(T, A)
    end
end

coord(b::BilinearMap3{N, U}) where {N, U <: Integer} = b.coord

base_ring(b::BilinearMap3{N, U}) where {N, U <: Integer} = b.base_ring

show(io::IO, b::BilinearMap3{N, U}) where {N, U <: Integer} = show(io, coord(b))

getindex(b::BilinearMap3{N, U}, i::Int) where {N, U <: Integer} = getindex(coord(b), i)

zero(b::BilinearMap3{N, U}) where {N, U <: Integer} = BilinearMap3(base_ring(b))

#=
function check_rings(a::BilinearMap3{N}, b::BilinearMap3{N}) where N
    base_ring(a) != base_ring(b) && error("Bilinear maps do not have the same base ring")
end
=#

function ==(a::BilinearMap3{N, U}, b::BilinearMap3{N, U}) where {N, U <: Integer}
#    check_rings(a, b)
    for j in 1:N
        if a[j] != b[j]
            return false
        end
    end
    return true
end

function ==(a::BilinearMap3{N, U}, i::Int) where {N, U <: Integer}
    i != 0 && error("Cannot compare a bilinear map with a nonzero integer")
    return iszero(a)
end

function add!(b::BilinearMap3{N, U}, decomp::Tuple{Vararg{Tuple{U, U}}},
              j::Int, d::Tuple{Vararg{nmod_mat}}) where {N, U <: Integer}

    zero!(b[j])
    p::U = modulus(base_ring(b[1]))

    for (x, c) in decomp
        M = -Int(c)*d[x + 0x01]
        t = translate(x, j, N, p)
        for i in j+1:N
            Nemo.add!(b[i], b[i], t[i-j]*M)
        end
    end
end

function copy(b::BilinearMap3{N, U}) where {N, U <: Integer}
    return BilinearMap3(deepcopy(coord(b)), base_ring(b)) # copy sufficient?
end

##############################################################################
#
# Type tryouts
#
# Simpler, better?
#
##############################################################################

struct ComAlgebra
    multiplication_tensor::Nemo.nmod_mat
    characteristic::Int
    degree::Int

    function ComAlgebra(p::Int, M::Nemo.nmod_mat)
        N = nrows(M)
        @assert ncols(M) == div(N*(N+1), 2)
        return new(M, p, N)
    end

    function ComAlgebra(k::FqNmodFiniteField)
        N = degree(k)
        p::Int = characteristic(k)
        M = MatrixSpace(ResidueRing(ZZ, p), N, div(N*(N+1), 2))()
        x = gen(k)
        for i in 1:N
            for j in i:N
                z = x^(i+j-2)
                col = div(N*(N+1)-(N+1-i)*(N+2-i), 2) + (j-i) + 1
                for l in 1:N
                    M[l, col] = coeff(z, l-1)
                end
            end
        end
        return new(M, p, N)
    end
end

multiplication_tensor(A::ComAlgebra) = A.multiplication_tensor
characteristic(A::ComAlgebra) = A.characteristic
degree(A::ComAlgebra) = A.degree

(A::ComAlgebra)() = MatrixSpace(base_ring(multiplication_tensor(A)), degree(A), 1)()
