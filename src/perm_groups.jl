# included from CGT_UniHeidelberg_2022.jl

const CGT = CGT_UniHeidelberg_2022
export PermutationGroup, order, gens, unsafe_gens, stabilizer_chain

abstract type Group end
abstract type AbstractPermGroup{P<:CGT.AbstractPermutation} <: Group end

mutable struct PermutationGroup{P,T<:CGT.AbstractTransversal} <: AbstractPermGroup{P}
    gens::Vector{P}
    order::BigInt
    stab_chain::CGT.StabilizerChain{T,P}

    # Constructor where:
    # only gens are known
    PermutationGroup{P, T}(gens::AbstractVector{P}) where {P,T} =
        new{P,T}(gens)
    # gens and order are known
    PermutationGroup{P, T}(gens::AbstractVector{P}, order::Integer) where {P,T} =
        new{P,T}(gens, order)

    # everything is known
    function PermutationGroup(
        gens::AbstractVector{P},
        order::Integer,
        stab_chain::CGT.StabilizerChain{T,P},
        check = true
    ) where {P,T}
        if check
            # we could/should add some consistency checks here e.g.
            @assert CGT.order(stab_chain) == order
            @assert all(gens) do g
                _, r = CGT.sift(stab_chain, g)
                isone(r)
            end
        end
        return new{P,T}(gens, order, stab_chain)
    end
end

PermutationGroup(gens::AbstractVector, sc::CGT.StabilizerChain{T}) where {T} =
    PermutationGroup(gens, order(sc), sc)
PermutationGroup(gens::AbstractVector{P}) where P = PermutationGroup{P, Transversal{Int64, P, typeof(^)}}(gens)

"""
    unsafe_gens(G::Group)
An unsafe version of `gens(G)`, the returned value may _alias_ internal data structures of `G`.

In particular should the returned value leave its caller scope, the safe version `gens(G)` must be used.
"""
unsafe_gens(G::Group) = G.gens
gens(G::Group) = copy(unsafe_gens(G))
gens(G::Group, i::Integer) = gens(G)[i]

# """
#     order([I=BigInt,] G::Group)
# Return order of group `G` as an instance of `I`.
# By default a `BigInt` (i.e. arbitrary sized integer) is returned.
# """
order(G::Group) = order(BigInt, G) # group orders can get very big very quickly
# order(sc::CGT.StabilizerChain) = order(BigInt, sc)

# order(::Type{I}, sc::CGT.StabilizerChain) where {I} =
#     convert(I, mapreduce(length, *, CGT.transversals(sc), init = one(I)))

"""
Lazy computation of stabilizer chain and group order.
"""
_knows_order(G::PermutationGroup) = isdefined(G, :order)

function order(::Type{I}, G::PermutationGroup) where {I<:Integer}
    if !_knows_order(G)
        G.order = order(BigInt, stabilizer_chain(G))
    end
    return convert(I, G.order)
end

function stabilizer_chain(G::PermutationGroup{P,T}) where {P,T}
    if !isdefined(G, :stab_chain)
        G.stab_chain = if _knows_order(G)
            CGT.schreier_sims(T, gens(G), order(G))
        else
            CGT.schreier_sims(T, gens(G))
        end
    end
    return G.stab_chain
end

"""
Other group methods and iteration protocol.
"""
function Base.in(G::AbstractPermGroup, p::CGT.AbstractPermutation)
    _, r = CGT.sift(stabilizer_chain(G), p)
    return isone(r)
end

Base.one(G::AbstractPermGroup) = one(first(gens(G)))
Base.eltype(::Type{<:AbstractPermGroup{P}}) where P = P
Base.length(G::Group) =
    order(G) > typemax(Int) ? typemax(Int) : order(Int, G) # practical limit when iterating
