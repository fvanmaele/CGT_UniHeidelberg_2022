# included from src/CGT_UniHeidelberg_2022.jl

const CGT = CGT_UniHeidelberg_2022

export decompose, base_images, sends_to

""" Express a group element `g` as a product of transversal elements, each of
    which has been obtained as a product of generators.
"""
function decompose(gens::AbstractVector{P}, g::P) where {P}
    # Use `SchreierTree` to keep track of how transversal elements arose as 
    # products of the original generators.
    G = PermutationGroup{P, CGT.SchreierTree{Int64, P, typeof(^)}}(gens)
    sc = CGT.stabilizer_chain(G)

    # Modified element check which keeps track of transversal elements.
    r = g # the residual
    L = Vector{P}()

    for i in 1:depth(sc)
        T = CGT.transversal(sc, i)
        y = CGT.basis(sc, i)^r

        if y ∉ T
            throw(ArgumentError("element is not in group"))
        else
            Ty, Ty_word = CGT.trace(T, y)
            r = r*inv(Ty)
            for gen in Ty_word
                push!(L, gen)
            end
            
            if r == one(r)
                return L
            end
        end
    end
    throw(ArgumentError("element is not in group"))
end

""" Modified element test which returns a list of base images. This function
    is inverse to `perm_by_images`.
"""
function base_images(sc::CGT.StabilizerChain, g::Permutation)
    r = g # the residual
    L = Vector{Int}()

    for i in 1:depth(sc)
        T = CGT.transversal(sc, i)
        y = CGT.basis(sc, i)^r

        if y ∉ T
            throw(ArgumentError("element is not in group"))
        else
            r = r*inv(T[y])
            push!(L, y)
        end
    end

    if r != one(r)
        throw(ArgumentError("element is not in group"))
    end
    @assert length(L) == CGT.depth(sc)
    return L
end

""" Function using backtrack search which, given 2 group elements `g` and `h`,
    finds a group element that sends `g` to `h`. To this end, the base images
    representing the group elements are retrieved. The search tree is traversed
    in breadth-first order, looking for a common ancestor of `g` and `h`.
"""
function sends_to(sc::CGT.StabilizerChain, g::CGT.Permutation, h::CGT.Permutation)

end