# included from src/CGT_UniHeidelberg_2022.jl

export perm_by_images
const CGT = CGT_UniHeidelberg_2022

""" Implement a method to reconstruct a group element from the images of a given
    basis of a `StabilizerChain`. This function is inverse to `base_images`.
"""
function perm_by_images(sc::StabilizerChain, γ::AbstractVector)
    @assert length(γ) == CGT.depth(sc)
    @assert length(γ) >= 1
    g = one(first(CGT.gens(sc, 1)))

    # proceed in opposite order of element test
    for i in CGT.depth(sc):-1:1
        T = CGT.transversal(sc, i)
        if γ[i] ∈ T
            g = g * T[γ[i]]
        else
            throw("group element not generated from base image")
        end
    end
    return g
end


""" Obtain uniformly random elements from the group generated by a `StabilizerChain`
"""
function Base.rand(sc::StabilizerChain)
    y = Int64[]
    for i in 1:CGT.depth(sc)
        T = CGT.transversal(sc, i)
        push!(y, rand(T)) # push random base image
    end
    g = perm_by_images(sc, y)
    return g
end
