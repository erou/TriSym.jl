#######################################################
#
# Search functions
#
#######################################################

export multiplication_search

function tri_symmetric_search_coord(b::Nemo.nmod_mat, E::AffineFieldElements,
         d::Dict{Nemo.fq_nmod, Nemo.nmod_mat}, L::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod,
         Int}}}, 1}, count::Int, bound::Int, G::Array{Tuple{Nemo.fq_nmod, Int},
         1} = Tuple{Nemo.fq_nmod, Int}[]) where N

    #println("coord: begin")
    #println(G)
    r = rank(b)
    p::Int = Nemo.base_ring(parent(b)).n

    if r <= bound
        if r == 0
            push!(L, Tuple(G))

        elseif count == 0
            if bound == 1
                println("case: bound = 1")
            end
            elems = Tuple{Nemo.fq_nmod, Int}[]
            for x in E
                for c in 1:p-1
                    if rank(b-c*d[x]) < r
                        push!(elems, (deepcopy(x), c))
                        break
                    end
                end
            end

            elems0 = Nemo.fq_nmod[x for (x, c) in elems]
            len = length(elems)
            for j in 1:len
                y, c = elems[j]
                B2 = elems0[j+1:end]
                G2 = copy(G)
                push!(G2, (y, c))
                tri_symmetric_search_coord(b-c*d[y], B2, d, L, bound-1, G2)
            end

        else
            for x in E
                y = next(E, x)
                E2 = AffineFieldElements(k, j, y)
                for c in 1:p-1
                    t2 = t-c*d[x]
                    r2 = rank(t2)
                    G2 = copy(G)
                    push!(G2, (deepcopy(x), c))
                    if r2 < r
                        tri_symmetric_search_coord(t2, E2, d, L, count, bound-1, G2)
                    elseif r2 == r
                        tri_symmetric_search_coord(t2, E2, d, L, count-1, bound-1, G2)
                    end
                end
            end
        end
    end
end

function tri_symmetric_search_coord(b::Nemo.nmod_mat, E::Array{fq_nmod, 1},
         d::Dict{Nemo.fq_nmod, Nemo.nmod_mat}, L::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod,
         Int}}}, 1}, bound::Int, G::Array{Tuple{Nemo.fq_nmod, Int}, 1}) where N

    #println("coord2: begin")
    #println(G)
    r = rank(b)
    p::Int = Nemo.base_ring(parent(b)).n

    if r <= bound
        if r == 0
            push!(L, Tuple(G))
        else
            elems = Tuple{Nemo.fq_nmod, Int}[]
            for x in E
                for c in 1:p-1
                    if rank(b-c*d[x]) < r
                        push!(elems, (deepcopy(x), c))
                        break
                    end
                end
            end

            elems0 = Nemo.fq_nmod[x for (x, c) in elems]
            len = length(elems)
            for j in 1:len
                y, c = elems[j]
                B2 = elems0[j+1:end]
                G2 = copy(G)
                push!(G2, (y, c))
                tri_symmetric_search_coord(b-c*d[y], B2, d, L, bound-1, G2)
            end
       end
    end
end

function tri_symmetric_search(B::BilinearMap{N}, d::Dict{Nemo.fq_nmod,
         Nemo.nmod_mat}, Lglob::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}},
         1}, counts::NTuple{N, Int}, bound::Int, j::Int = 1,
         Lloc::Array{Tuple{Nemo.fq_nmod, Int}, 1} = Tuple{Nemo.fq_nmod, Int}[]) where N

    println("B, Lglob, j = ")
    println(B)
    println(Lglob)
    println(j)
    if iszero(B)
        push!(Lglob, Tuple(Lloc))
    elseif rank(B[j]) <= bound
        M = Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}[]
        count = counts[j]

        E = AffineFieldElements(base_ring(B), j)
        tri_symmetric_search_coord(B[j], E, d, M, count, bound)

        for m in M
            B2 = deepcopy(B)
            add!(B2, m, j, d)
            Lloc2 = copy(Lloc)
            append!(Lloc2, m)
            tri_symmetric_search(B2, d, Lglob, counts, bound-length(m), j+1, Lloc2)
        end
    end
end

function multiplication_search(k::Nemo.FqNmodFiniteField, counts::NTuple{N,
                               Int}, bound::Int) where N
    B = multiplication_bilinear_map(k)
    d = make_conversion_dict(k)
    Lglob = Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}[]
    tri_symmetric_search(B, d, Lglob, counts, bound)
    return Lglob
end 
