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
        println("depth: ", depth, " base image: ", δ)
        if CGT.depth(sc) == depth  # we are in a leaf node
            push!(L, T[δ]*g)
        else
            backtrack!(sc, L; g=T[δ]*g, depth=depth+1)
        end
    end
    return L
end

""" Version using `Channel` instead of a vector, for use with the
    iterator interface.
"""
function backtrack!(sc::CGT.StabilizerChain, C::Channel;
                    g::CGT.Permutation = CGT.Permutation(Int[]), depth::Int = 1)
    T = CGT.transversal(sc, depth)

    for δ ∈ T
        if CGT.depth(sc) == depth  # we are in a leaf node
            put!(C, T[δ]*g)
        else
            backtrack!(sc, C; g=T[δ]*g, depth=depth+1)
        end
    end
end

""" Experimental version using explicit stack.
    Note: the order is breadth-first, not depth-first. This can be used
    to easily prune the search tree using partial base images.
"""
function backtrack_stack!(sc::CGT.StabilizerChain, L::AbstractVector)
    stack = [(1, CGT.Permutation(Int[]))]

    while !isempty(stack)
        depth, g = pop!(stack)
        T = CGT.transversal(sc, depth)

        for δ ∈ T
            #println("depth: ", depth, " base image: ", δ)
            if CGT.depth(sc) == depth  # we are in a leaf node
                #println("group element: ", g*T[δ])
                push!(L, T[δ]*g)
            else
                push!(stack, (depth+1, T[δ]*g))
            end
        end
    end
    return L
end

""" Helpers to make backtrack! work with the iterator interface.
"""
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

# XXX: due to a race condition with channels, Base.collect returns a vector
# with undefined elements.
function Base.collect(G::PermutationGroup)
    n = order(G)
    L = Vector{eltype(G)}(undef, n)
    state = PGroupIterator(G)
    for i = 1:n
        g = take!(state)
        L[i] = g
    end
    return L
end
