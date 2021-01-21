#############################################################################
#
# Search functions
# ================
#
# Algorithms searching for trisymmetric decompositions using the vector space
# structure of the algebra. The algorithms are essentially all the same, with
# minor differences in how the data is represented, resulting in different
# behaviour/speed in practice.
#
#############################################################################

export multiplication_search

#############################################################################
#
# First search function
#
#############################################################################

function tri_symmetric_search_coord(b::Nemo.nmod_mat, E::Array{Nemo.fq_nmod, 1},
         d::Dict{Nemo.fq_nmod, Nemo.nmod_mat}, L::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod,
         Int}}}, 1}, count::Int, bound::Int, G::Array{Tuple{Nemo.fq_nmod, Int},
         1} = Tuple{Nemo.fq_nmod, Int}[]) where N

#    println("G, count, bound, E =")
#    println(G)
#    println(count)
#    println(bound)
#    println(E)
    r = rank(b)
    p::UInt = Nemo.base_ring(parent(b)).n

    if r <= bound
        if r == 0
            push!(L, Tuple(G))

        elseif count == 0
            b2 = similar(b)
            elems = Tuple{Nemo.fq_nmod, Int, Nemo.nmod_mat}[]
            for x in E
                for c in 1:p-1
                #    #if bound == 1
                #        println("case: bound = 1")
                #        return (d, x)
                #    end
                    submul!(b2, b, d[x], c)
                    if rank(b2) < r
                        push!(elems, (x, c, Base.copy(b2)))
                        break
                    end
                end
            end

            elems0 = Nemo.fq_nmod[x[1] for x in elems]
            for j in 1:length(elems)
                y, c, b2 = elems[j]
                B2 = elems0[j+1:end]
                G2 = Base.copy(G) #not **deep** copy, be careful!
                push!(G2, (y, c))
                tri_symmetric_search_coord(b2, B2, d, L, bound-1, G2)
            end

        else
            b2 = similar(b)
            for j in 1:length(E)
                x = E[j]
                E2 = E[j+1:end]
                for c in 1:p-1
                    submul!(b2, b, d[x], c) # no alloc for b2 here
                    r2 = rank(b2)
                    G2 = Base.copy(G) #not **deep** copy, be careful!
                    push!(G2, (x, c))
                    if r2 < r
                        tri_symmetric_search_coord(b2, E2, d, L, count, bound-1, G2)
                    elseif r2 == r
                        tri_symmetric_search_coord(b2, E2, d, L, count-1, bound-1, G2)
                    end
                end
            end
        end
    end
end

function tri_symmetric_search_coord(b::Nemo.nmod_mat, E::Array{fq_nmod, 1},
         d::Dict{Nemo.fq_nmod, Nemo.nmod_mat}, L::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod,
         Int}}}, 1}, bound::Int, G::Array{Tuple{Nemo.fq_nmod, Int}, 1}) where N

#    println("* G, bound, E =")
#    println(G)
#    println(bound)
#    println(E)
    r = rank(b)
    p::UInt = Nemo.base_ring(parent(b)).n

    if r <= bound
        if r == 0
            push!(L, Tuple(G))
        else
            b2 = similar(b)
            elems = Tuple{Nemo.fq_nmod, Int, Nemo.nmod_mat}[]
            for x in E
                for c in 1:p-1
                    submul!(b2, b, d[x], c)
                    if rank(b2) < r
                        push!(elems, (x, c, Base.copy(b2)))
                        break
                    end
                end
            end

            elems0 = Nemo.fq_nmod[x[1] for x in elems]
            for j in 1:length(elems)
                y, c, b2 = elems[j]
                B2 = elems0[j+1:end]
                G2 = Base.copy(G) #not **deep** copy, be careful!
                push!(G2, (y, c))
                tri_symmetric_search_coord(b2, B2, d, L, bound-1, G2)
            end
       end
    end
end

function tri_symmetric_search(B::BilinearMap{N}, d::Dict{Nemo.fq_nmod,
         Nemo.nmod_mat}, Lglob::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}},
         1}, counts::NTuple{N, Int}, bound::Int, j::Int = 1,
         Lloc::Array{Tuple{Nemo.fq_nmod, Int}, 1} = Tuple{Nemo.fq_nmod, Int}[]) where N

#    println("*****")
#    println("B, Lloc, j = ")
#    println(B)
#    println(Lloc)
#    println(j)
#    println("*****")
    if iszero(B)
        push!(Lglob, Tuple(Lloc))
    elseif rank(B[j]) <= bound
        M = Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}[]
        count = counts[j]

        E = AffineFieldElements2(base_ring(B), j)
        tri_symmetric_search_coord(B[j], E, d, M, count, bound)
        #res = tri_symmetric_search_coord(B[j], E, d, M, count, bound)
        #if res != nothing
        #    println("returned")
        #    return res
        #end

        for m in M
            B2 = copy(B)
            add!(B2, m, j, d)
            Lloc2 = Base.copy(Lloc) #not **deep** copy, be careful!
            append!(Lloc2, m)
            tri_symmetric_search(B2, d, Lglob, counts, bound-length(m), j+1, Lloc2)
            #res = tri_symmetric_search(B2, d, Lglob, counts, bound-length(m), j+1, Lloc2)
            #if res != nothing
            #    println("returned -")
            #    return res
            #end
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

#############################################################################
#
# Search function with time/memory tradeoff
#
#############################################################################

function tri_symmetric_search_coord_memo(b::Nemo.nmod_mat, E::Array{Nemo.fq_nmod, 1},
         d::Dict{Nemo.fq_nmod, Nemo.nmod_mat}, L::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod,
         Int}}}, 1}, count::Int, bound::Int, precalc::Dict{Nemo.nmod_mat,
         Tuple{Vararg{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}}}}, precalc_bound::Int,
         G::Array{Tuple{Nemo.fq_nmod, Int}, 1} = Tuple{Nemo.fq_nmod, Int}[]) where N

    #println("G, count, bound, E =")
    #println(G)
    #println(count)
    #println(bound)
    #println(E)
    r = rank(b)
    p::UInt = Nemo.base_ring(parent(b)).n

    if r <= bound
        if r <= precalc_bound
            if haskey(precalc, b)
                decomps = precalc[b]
                for dec in decomps
                    if length(dec) <= bound
                        G2 = Base.copy(G)
                        for y in dec
                            push!(G2, y)
                        end
                        sort!(G2)
                        T2 = Tuple(G2)
                        if !(T2 in L)
                            push!(L, T2)
                        end
                    end
                end
            end
            if count > 0
                b2 = similar(b)
                for j in 1:length(E)
                    x = E[j]
                    E2 = E[j+1:end]
                    for c in 1:p-1
                        submul!(b2, b, d[x], c) # no alloc for b2 here
                        r2 = rank(b2)
                        G2 = Base.copy(G) #not **deep** copy, be careful!
                        push!(G2, (x, c))
                        if r2 == r
                            tri_symmetric_search_coord_memo(b2, E2, d, L, count-1,
                                                       bound-1, precalc, precalc_bound, G2)
                        end
                    end
                end
            end

        elseif count == 0
            b2 = similar(b)
            elems = Tuple{Nemo.fq_nmod, Int, Nemo.nmod_mat}[]
            for x in E
                for c in 1:p-1
                #    #if bound == 1
                #        println("case: bound = 1")
                #        return (d, x)
                #    end
                    submul!(b2, b, d[x], c)
                    if rank(b2) < r
                        push!(elems, (x, c, Base.copy(b2)))
                        break
                    end
                end
            end

            elems0 = Nemo.fq_nmod[x[1] for x in elems]
            for j in 1:length(elems)
                y, c, b2 = elems[j]
                B2 = elems0[j+1:end]
                G2 = Base.copy(G) #not **deep** copy, be careful!
                push!(G2, (y, c))
                tri_symmetric_search_coord_memo(b2, B2, d, L, bound-1, precalc,
                                                precalc_bound, G2)
            end

        else
            b2 = similar(b)
            for j in 1:length(E)
                x = E[j]
                E2 = E[j+1:end]
                for c in 1:p-1
                    submul!(b2, b, d[x], c) # no alloc for b2 here
                    r2 = rank(b2)
                    G2 = Base.copy(G) #not **deep** copy, be careful!
                    push!(G2, (x, c))
                    if r2 < r
                        tri_symmetric_search_coord_memo(b2, E2, d, L, count, bound-1,
                                                        precalc, precalc_bound, G2)
                    elseif r2 == r
                        tri_symmetric_search_coord_memo(b2, E2, d, L, count-1, bound-1,
                                                        precalc, precalc_bound, G2)
                    end
                end
            end
        end
    end
end

function tri_symmetric_search_coord_memo(b::Nemo.nmod_mat, E::Array{fq_nmod, 1},
         d::Dict{Nemo.fq_nmod, Nemo.nmod_mat}, L::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod,
         Int}}}, 1}, bound::Int, precalc::Dict{Nemo.nmod_mat,
         Tuple{Vararg{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}}}},
         precalc_bound::Int, G::Array{Tuple{Nemo.fq_nmod, Int}, 1}) where N

    #println("* G, bound, E =")
    #println(G)
    #println(bound)
    #println(E)

    r = rank(b)
    p::UInt = Nemo.base_ring(parent(b)).n

    if r <= bound
        if r <= precalc_bound
            if haskey(precalc, b)
                decomps = precalc[b]
                for dec in decomps
                    if length(dec) <= bound
                        G2 = Base.copy(G)
                        for y in dec
                            push!(G2, y)
                        end
                        sort!(G2)
                        T2 = Tuple(G2)
                        if !(T2 in L)
                            push!(L, T2)
                        end
                    end
                end
            end

        else
            b2 = similar(b)
            elems = Tuple{Nemo.fq_nmod, Int, Nemo.nmod_mat}[]
            for x in E
                for c in 1:p-1
                #    #if bound == 1
                #        println("case: bound = 1")
                #        return (d, x)
                #    end
                    submul!(b2, b, d[x], c)
                    if rank(b2) < r
                        push!(elems, (x, c, Base.copy(b2)))
                        break
                    end
                end
            end

            elems0 = Nemo.fq_nmod[x[1] for x in elems]
            for j in 1:length(elems)
                y, c, b2 = elems[j]
                B2 = elems0[j+1:end]
                G2 = Base.copy(G) #not **deep** copy, be careful!
                push!(G2, (y, c))
                tri_symmetric_search_coord_memo(b2, B2, d, L, bound-1, precalc,
                                                precalc_bound, G2)
            end
        end
    end
end

function tri_symmetric_search_memo(B::BilinearMap{N}, d::Dict{Nemo.fq_nmod,
         Nemo.nmod_mat}, Lglob::Array{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}},
         1}, counts::NTuple{N, Int}, bound::Int, precalc::Array{Dict{Nemo.nmod_mat,
         Tuple{Vararg{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}}}}, 1},
         precalc_bound::Int, j::Int = 1,
         Lloc::Array{Tuple{Nemo.fq_nmod, Int}, 1} = Tuple{Nemo.fq_nmod, Int}[]) where N

    #println("*****")
    #println("B, Lloc, j = ")
    #println(B)
    #println(Lloc)
    #println(j)
    #println("*****")
    if iszero(B)
        push!(Lglob, Tuple(Lloc))
    elseif rank(B[j]) <= bound
        M = Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}[]
        count = counts[j]
        precalc_dict = precalc[j]

        E = AffineFieldElements2(base_ring(B), j)
        tri_symmetric_search_coord_memo(B[j], E, d, M, count, bound,
                                        precalc_dict, precalc_bound)
        #res = tri_symmetric_search_coord(B[j], E, d, M, count, bound)
        #if res != nothing
        #    println("returned")
        #    return res
        #end

        for m in M
            B2 = copy(B)
            add!(B2, m, j, d)
            Lloc2 = Base.copy(Lloc) #not **deep** copy, be careful!
            append!(Lloc2, m)
            tri_symmetric_search_memo(B2, d, Lglob, counts, bound-length(m),
                                      precalc, precalc_bound, j+1, Lloc2)
            #res = tri_symmetric_search(B2, d, Lglob, counts, bound-length(m), j+1, Lloc2)
            #if res != nothing
            #    println("returned -")
            #    return res
            #end
        end
    end
end

function multiplication_search_memo(k::Nemo.FqNmodFiniteField, counts::NTuple{N,
                               Int}, bound::Int, precalc_bound::Int) where N
    B = multiplication_bilinear_map(k)
    d = make_conversion_dict(k)
    precalc = Array{Dict{Nemo.nmod_mat,
               Tuple{Vararg{Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}}}}, 1}(undef, N)

    for i in 1:N
        precalc[i] = make_small_trisym_coord(k, precalc_bound, d, i)
    end

    Lglob = Tuple{Vararg{Tuple{Nemo.fq_nmod, Int}}}[]
    
    tri_symmetric_search_memo(B, d, Lglob, counts, bound, precalc, precalc_bound)
    return Lglob
end 

#############################################################################
#
# Search function for AbstractCommutativeAlgebra
#
#############################################################################

function tri_symmetric_search_coord(b::nmod_mat, E::Array{nmod_mat, 1},
         d::Dict{nmod_mat, nmod_mat}, L::Array{Tuple{Vararg{Tuple{nmod_mat,
         Int}}}, 1}, count::Int, bound::Int, G::Array{Tuple{nmod_mat, Int},
         1} = Tuple{nmod_mat, Int}[]) where N

    #println("G, count, bound, E =")
    #println(G)
    #println(count)
    #println(bound)
    #println(E)
    r = rank(b)
    p::UInt = Nemo.base_ring(parent(b)).n

    if r <= bound
        if r == 0
            push!(L, Tuple(G))

        elseif count == 0
            b2 = similar(b)
            elems = Tuple{nmod_mat, Int, nmod_mat}[]
            for x in E
                for c in 1:p-1
                #    #if bound == 1
                #        println("case: bound = 1")
                #        return (d, x)
                #    end
                    submul!(b2, b, d[x], c)
                    if rank(b2) < r
                        push!(elems, (x, c, copy(b2)))
                        break
                    end
                end
            end

            elems0 = nmod_mat[x[1] for x in elems]
            for j in 1:length(elems)
                y, c, b2 = elems[j]
                B2 = elems0[j+1:end]
                G2 = Base.copy(G) #not **deep** copy, be careful!
                push!(G2, (y, c))
                tri_symmetric_search_coord(b2, B2, d, L, bound-1, G2)
            end

        else
            b2 = similar(b)
            for j in 1:length(E)
                x = E[j]
                E2 = E[j+1:end]
                for c in 1:p-1
                    submul!(b2, b, d[x], c) # no alloc for b2 here
                    r2 = rank(b2)
                    G2 = Base.copy(G) #not **deep** copy, be careful!
                    push!(G2, (x, c))
                    if r2 < r
                        tri_symmetric_search_coord(b2, E2, d, L, count, bound-1, G2)
                    elseif r2 == r
                        tri_symmetric_search_coord(b2, E2, d, L, count-1, bound-1, G2)
                    end
                end
            end
        end
    end
end

function tri_symmetric_search_coord(b::nmod_mat, E::Array{nmod_mat, 1},
         d::Dict{nmod_mat, nmod_mat}, L::Array{Tuple{Vararg{Tuple{nmod_mat,
         Int}}}, 1}, bound::Int, G::Array{Tuple{nmod_mat, Int}, 1}) where N

    #println("* G, bound, E =")
    #println(G)
    #println(bound)
    #println(E)
    r = rank(b)
    p::UInt = Nemo.base_ring(parent(b)).n

    if r <= bound
        if r == 0
            push!(L, Tuple(G))
        else
            b2 = similar(b)
            elems = Tuple{nmod_mat, Int, Nemo.nmod_mat}[]
            for x in E
                for c in 1:p-1
                    submul!(b2, b, d[x], c)
                    if rank(b2) < r
                        push!(elems, (x, c, Base.copy(b2)))
                        break
                    end
                end
            end

            elems0 = Nemo.nmod_mat[x[1] for x in elems]
            for j in 1:length(elems)
                y, c, b2 = elems[j]
                B2 = elems0[j+1:end]
                G2 = Base.copy(G) #not **deep** copy, be careful!
                push!(G2, (y, c))
                tri_symmetric_search_coord(b2, B2, d, L, bound-1, G2)
            end
       end
    end
end

function tri_symmetric_search(B::BilinearMap2{N}, d::Dict{nmod_mat,
         nmod_mat}, Lglob::Array{Tuple{Vararg{Tuple{nmod_mat, Int}}},
         1}, counts::NTuple{N, Int}, bound::Int, j::Int = 1,
         Lloc::Array{Tuple{nmod_mat, Int}, 1} = Tuple{nmod_mat, Int}[]) where N

    #println("*****")
    #println("B, Lloc, j = ")
    #println(B)
    #println(Lloc)
    #println(j)
    #println("*****")
    if iszero(B)
        push!(Lglob, Tuple(Lloc))
    elseif rank(B[j]) <= bound
        M = Tuple{Vararg{Tuple{nmod_mat, Int}}}[]
        count = counts[j]

        E = AffineFieldElements2(base_ring(B), j)
        tri_symmetric_search_coord(B[j], E, d, M, count, bound)
        #res = tri_symmetric_search_coord(B[j], E, d, M, count, bound)
        #if res != nothing
        #    println("returned")
        #    return res
        #end

        for m in M
            B2 = copy(B)
            add!(B2, m, j, d)
            Lloc2 = Base.copy(Lloc) #not **deep** copy, be careful!
            append!(Lloc2, m)
            tri_symmetric_search(B2, d, Lglob, counts, bound-length(m), j+1, Lloc2)
            #res = tri_symmetric_search(B2, d, Lglob, counts, bound-length(m), j+1, Lloc2)
            #if res != nothing
            #    println("returned -")
            #    return res
            #end
        end
    end
end

function multiplication_search(A::AbstractCommutativeAlgebra, L::nmod_mat, counts::NTuple{N,
                               Int}, bound::Int) where N
    B = multiplication_bilinear_map(A)
    d = make_conversion_dict(A, L)
    Lglob = Tuple{Vararg{Tuple{nmod_mat, Int}}}[]
    
    tri_symmetric_search(B, d, Lglob, counts, bound)
    return Lglob
end 


#############################################################################
#
# Search function for CommutativeAlgebra
#
#############################################################################

function tri_symmetric_search_coord(b::nmod_mat, E::Tuple{Vararg{U}},
                                    d::Tuple{Vararg{nmod_mat}},
                                    L::Array{Tuple{Vararg{Tuple{U,
         U}}}, 1}, count::Int, bound::Int, G::Array{Tuple{U, U},
         1} = Tuple{U, U}[]) where {N, U <: Integer}

    #println("G, count, bound, E =")
    #println(G)
    #println(count)
    #println(bound)
    #println(E)
    r = rank(b)
    p::U = modulus(base_ring(parent(b)))
    o = one(U)

    if r <= bound
        if r == 0
            push!(L, Tuple(G))

        elseif count == 0
            b2 = similar(b)
            elems = Tuple{U, U, nmod_mat}[]
            for x in E
                for c in o:p-1
                #    #if bound == 1
                #        println("case: bound = 1")
                #        return (d, x)
                #    end
                    submul!(b2, b, d[x + 0x01], UInt(c))
                    if rank(b2) < r
                        push!(elems, (x, c, copy(b2)))
                        break
                    end
                end
            end

            elems0 = Tuple(U[x[1] for x in elems])
            for j in 1:length(elems)
                y, c, b2 = elems[j]
                B2 = elems0[j+1:end]
                G2 = Base.copy(G) #not **deep** copy, be careful!
                push!(G2, (y, c))
                tri_symmetric_search_coord(b2, B2, d, L, bound-1, G2)
            end

        else
            b2 = similar(b)
            for j in 1:length(E)
                x = E[j]
                E2 = E[j+1:end]
                for c in o:p-1
                    submul!(b2, b, d[x + 0x01], UInt(c)) # no alloc for b2 here
                    r2 = rank(b2)
                    G2 = Base.copy(G) #not **deep** copy, be careful!
                    push!(G2, (x, c))
                    if r2 < r
                        tri_symmetric_search_coord(b2, E2, d, L, count, bound-1, G2)
                    elseif r2 == r
                        tri_symmetric_search_coord(b2, E2, d, L, count-1, bound-1, G2)
                    end
                end
            end
        end
    end
end

function tri_symmetric_search_coord(b::nmod_mat, E::Tuple{Vararg{U}},
         d::Tuple{Vararg{nmod_mat}}, L::Array{Tuple{Vararg{Tuple{U,
         U}}}, 1}, bound::Int, G::Array{Tuple{U, U}, 1}) where {N, U <: Integer}

    #println("* G, bound, E =")
    #println(G)
    #println(bound)
    #println(E)
    r = rank(b)
    p::U = modulus(base_ring(parent(b)))
    o = one(U)

    if r <= bound
        if r == 0
            push!(L, Tuple(G))
        else
            b2 = similar(b)
            elems = Tuple{U, U, Nemo.nmod_mat}[]
            for x in E
                for c in o:p-1
                    submul!(b2, b, d[x + 0x01], UInt(c))
                    if rank(b2) < r
                        push!(elems, (x, c, Base.copy(b2)))
                        break
                    end
                end
            end

            elems0 = Tuple(U[x[1] for x in elems])
            for j in 1:length(elems)
                y, c, b2 = elems[j]
                B2 = elems0[j+1:end]
                G2 = Base.copy(G) #not **deep** copy, be careful!
                push!(G2, (y, c))
                tri_symmetric_search_coord(b2, B2, d, L, bound-1, G2)
            end
       end
    end
end

function tri_symmetric_search(B::BilinearMap3{N}, d::Tuple{Vararg{nmod_mat}},
         Lglob::Array{Tuple{Vararg{Tuple{U, U}}}, 1},
         counts::NTuple{N, Int}, bound::Int, j::Int = 1,
         Lloc::Array{Tuple{U, U}, 1} = Tuple{U, U}[]) where {N, U <: Integer}

    #println("*****")
    #println("B, Lloc, j = ")
    #println(B)
    #println(Lloc)
    #println(j)
    #println("*****")
    if iszero(B)
        push!(Lglob, Tuple(Lloc))
    elseif rank(B[j]) <= bound
        M = Tuple{Vararg{Tuple{U, U}}}[]
        count = counts[j]

        E = AffineElements(base_ring(B), U(j))
        tri_symmetric_search_coord(B[j], E, d, M, count, bound)
        #res = tri_symmetric_search_coord(B[j], E, d, M, count, bound)
        #if res != nothing
        #    println("returned")
        #    return res
        #end

        for m in M
            B2 = copy(B)
            add!(B2, m, j, d)
            Lloc2 = Base.copy(Lloc) #not **deep** copy, be careful!
            append!(Lloc2, m)
            tri_symmetric_search(B2, d, Lglob, counts, bound-length(m), j+1, Lloc2)
            #res = tri_symmetric_search(B2, d, Lglob, counts, bound-length(m), j+1, Lloc2)
            #if res != nothing
            #    println("returned -")
            #    return res
            #end
        end
    end
end
#=
function multiplication_search(A::AbstractCommutativeAlgebra, L::nmod_mat, counts::NTuple{N,
                               Int}, bound::Int) where N
    B = multiplication_bilinear_map(A)
    d = make_conversion_dict(A, L)
    Lglob = Tuple{Vararg{Tuple{nmod_mat, Int}}}[]
    
    tri_symmetric_search(B, d, Lglob, counts, bound)
    return Lglob
end 
=#

#############################################################################
#
# Search functions with vector space strategy
#
#############################################################################

function tri_symmetric_vector_search(W::nmod_mat, G::nmod_mat, d::Dict{nmod_mat,
                                     nmod_mat}, L::Array{nmod_mat, 1},
                                     elems::Array{nmod_mat, 1}, bound::Int,
                                     count::Int)

    if count < bound
        for (j, x) in enumerate(elems)
            W2 = copy(W)
            add_vector!(W2, x)
            elems2 = mycopy(elems[j+1:end])
            reduce_column!(elems2, x)
            tri_symmetric_vector_search(W2, G, d, L, elems2, bound, count+1)
        end
    elseif count == bound
        if rank(intersection(W, G)) == bound
            push!(L, W)
        end
    end
end

function multiplication_vector_search(A::AbstractCommutativeAlgebra{N},
                                      L::nmod_mat, bound::Int) where N
    d = make_conversion_dict_vectors(A, L)
    W = multiplication_vector(A)
    G = construct_matrix(d, A)
    M = nmod_mat[]

    elems = Array{nmod_mat, 1}(undef, ncols(G))
    for j in 1:ncols(G)
        elems[j] = G[:, j]
    end

    return W, G, d, M, elems, bound

    tri_symmetric_vector_search(W, G, d, M, elems, bound, 1)

    return M
end
