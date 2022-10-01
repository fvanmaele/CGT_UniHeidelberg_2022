# included from src/CGT_UniHeidelberg_2022.jl

const CGT = CGT_UniHeidelberg_2022

export split_to_gens

# Express a group element `g` as a product of transversal elements, each of
# which has been obtained as a product of generators.
function split_to_gens(gens::AbstractVector{P}, g::P) where {P}
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
