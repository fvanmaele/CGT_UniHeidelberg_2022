# included from CGT_UniHeidelberg_2022.jl

const CGT = CGT_UniHeidelberg_2022

export backtrack!, backtrack_stack!

""" Basic version of backtrack search without oracle.
    Returns a list of all elements of G.
"""
function backtrack!(L::AbstractVector, sc::CGT.StabilizerChain;
                    g::CGT.Permutation = CGT.Permutation(Int[]), depth::Int = 1)
    T = CGT.transversal(sc, depth)

    for δ ∈ T
        if CGT.depth(sc) == depth  # we are in a leaf node
            push!(L, g*T[δ])
        else
            backtrack!(L, sc, g=g*T[δ], depth=depth+1)
        end
    end
    return L
end

# XXX: Version using stack for experimentation purposes
function backtrack_stack!(L::AbstractVector, sc::CGT.StabilizerChain)
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
