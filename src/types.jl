#######################################################
#
# Type to represent bilinear maps
#
#######################################################

export BilinearMap

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

Base.show(io::IO, b::BilinearMap{N}) where N = show(io, coord(b))

Base.getindex(b::BilinearMap{N}, i::Int) where N = getindex(coord(b), i)

Base.setindex!(b::BilinearMap{N}, i::Int) where N = getindex(coord(b), i)

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
