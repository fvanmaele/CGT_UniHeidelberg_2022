export schreier_sims, StabilizerChain

struct StabilizerChain{T<:AbstractTransversal, P<:AbstractPermutation}
    transversals::Vector{T}
    gens::Vector{Vector{P}} # generators for stabilizers

    function StabilizerChain(
        transversals::AbstractVector{<:AbstractTransversal},
        gens::AbstractVector{<:AbstractVector{<:AbstractPermutation}},
        # for some reason I could only get the tests to work when I specified the
        # type here, so using 'Permutation' instead of '<:AbstractPermutation'
    )
        @assert length(transversals) == length(gens)
        T = eltype(transversals)
        P = eltype(eltype(gens))
        return new{T,P}(transversals, gens)
    end
end

StabilizerChain{T, P}() where {T,P} = StabilizerChain(Vector{T}(), Vector{Vector{P}}())

depth(sc::StabilizerChain) = length(sc.transversals)
transversals(sc::StabilizerChain) = sc.transversals
transversal(sc::StabilizerChain, i::Integer) = sc.transversals[i]
gens(sc::StabilizerChain, i::Integer) = sc.gens[i]

basis(sc::StabilizerChain, i::Integer) = first(transversal(sc, i))
basis(sc::StabilizerChain) = [basis(sc, i) for i in 1:depth(sc)]

schreier_sims(S::AbstractVector{<:AbstractPermutation}) =
    schreier_sims(Transversal{Int64, eltype(S), typeof(^)}, S)

function schreier_sims(::Type{T}, S::AbstractVector{P}) where {T<:AbstractTransversal, P<:AbstractPermutation}
    sc = StabilizerChain{T, P}()

    for s in S
        push!(sc, s)
    end

    return sc
end

function sift(sc::StabilizerChain, g::AbstractPermutation; start_at_depth=1)
    r = g # the residual

    for i in start_at_depth:depth(sc)
        T = transversal(sc, i)
        β = first(T)
        γ = action(T)(β, r)

        if γ ∉ T
            return i, r
        else
            r = r*inv(T[γ])
            if r == one(g)
                return i+1, r
            end
        end
    end
    return depth(sc)+1, r
end

function Base.push!(sc::StabilizerChain, g::AbstractPermutation; at_depth::Integer=1)
    isone(g) && return sc
    @assert 1 ≤ at_depth ≤ depth(sc)+1

    for i in 1:at_depth-1
        β = basis(sc, i)
        @assert action(transversal(sc,i))(β, g) == β
    end

    _, r = sift(sc, g, start_at_depth=at_depth)

    isone(r) && return sc # g is already in the group generated by sc

    if at_depth == depth(sc) + 1
        extend_chain!(sc, r)
    else
        extend_gens!(sc, r, at_depth=at_depth)
    end
    return sc
end

function extend_chain!(sc::StabilizerChain{T}, g::AbstractPermutation) where T
    β = AbstractPermutations.firstmoved(g)
    S = [g]
    tr = T(β, S, ^)
    push!(sc.transversals, tr)
    push!(sc.gens, S)
    # a potentially new Schreier generator on the next layer:
    schr = g^length(tr)
    # clearly β^schr == β, so it is a generator for the next stabilizer
    if !isone(schr)
        # here we can directly extend_chain and shortcut all checks + call to sift
        extend_chain!(sc, schr)
        # otherwise we could do
        # push!(sc, schr, at_depth=depth(sc)+1)
    end
    return sc
end

function extend_gens!(sc::StabilizerChain, g::AbstractPermutation; at_depth::Integer)
    push!(sc.gens[at_depth], g) # add new generators
    T = transversal(sc, at_depth) # in-place update of stabilizer chain

    # old points only with new generator
    for δ ∈ T
        γ = action(T)(δ, g)
        if γ ∉ T
            # add γ to orbit, update transversal
            push!(T, γ=>T[δ]*g)
        else
            # γ is a coset representative for G(i+1) in G(i)
            s = T[δ]*g*inv(T[γ]) # Schreier generator
            push!(sc, s, at_depth=at_depth+1) # push to next layer
        end
    end

    # recompute transversal + process all the new schreier generators
    for δ ∈ T
        for b ∈ sc.gens[at_depth] # includes g
            γ = action(T)(δ, b)
            if γ ∉ T
                push!(T, γ=>T[δ]*b)
            else
                s = T[δ]*b*inv(T[γ])
                push!(sc, s, at_depth=at_depth+1)
            end
        end
    end
    return sc
end

order(sc::StabilizerChain) = order(BigInt, sc)
order(::Type{I}, sc::StabilizerChain) where I =
    convert(I, mapreduce(length, *, transversals(sc), init=one(I)))
