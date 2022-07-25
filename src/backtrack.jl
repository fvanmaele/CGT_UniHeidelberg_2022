# included from CGT_UniHeidelberg_2022.jl

const CGT = CGT_UniHeidelberg_2022

export backtrack!, backtrack_stack!

""" Basic version of backtrack search without oracle.
    Returns a list of all elements of G.
"""
function backtrack!(sc::CGT.StabilizerChain, L::AbstractVector;
                    g::CGT.Permutation = CGT.Permutation(Int[]), depth::Int = 1)
    T = CGT.transversal(sc, depth)

    for δ ∈ T
        if CGT.depth(sc) == depth  # we are in a leaf node
            push!(L, g*T[δ])
        else
            backtrack!(sc, L, g=g*T[δ], depth=depth+1)
        end
    end
    return L
end

function backtrack!(sc::CGT.StabilizerChain, C::Channel;
                    g::CGT.Permutation = CGT.Permutation(Int[]), depth::Int = 1)
    T = CGT.transversal(sc, depth)

    for δ ∈ T
        if CGT.depth(sc) == depth  # we are in a leaf node
            put!(C, g*T[δ])
        else
            backtrack!(sc, C, g=g*T[δ], depth=depth+1)
        end
    end
end

# Experimental version using explicit stack
# XXX: the ordering is different from the recursive version
function backtrack_stack!(sc::CGT.StabilizerChain, L::AbstractVector)
    stack = [(1, CGT.Permutation(Int[]))]

    while !isempty(stack)
        depth, g = pop!(stack)
        T = CGT.transversal(sc, depth)

        for δ ∈ T
            if CGT.depth(sc) == depth  # we are in a leaf node
                push!(L, g*T[δ])
            else
                push!(stack, (depth+1, g*T[δ]))
            end
        end
    end
    return L
end

# Helpers to make backtrack! work with iterator interface
struct PGroupIterator
    C::Channel

    function PGroupIterator(G::CGT.PermutationGroup)
        C = Channel((channel_arg) ->
            backtrack!(CGT.stabilizer_chain(G), channel_arg))
        return new(C)
    end
end

Base.take!(It::PGroupIterator) = take!(It.C)
Base.isready(It::PGroupIterator) = isready(It.C)

function Base.iterate(G::CGT.PermutationGroup)
    state = PGroupIterator(G)
    return (take!(state), state) # should contain at least one element
end

function Base.iterate(::CGT.PermutationGroup, state::PGroupIterator)
    if isready(state)
        return (take!(state), state)
    else
        return nothing
    end
end
