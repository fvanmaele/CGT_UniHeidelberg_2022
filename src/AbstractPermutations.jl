module AbstractPermutations

import ..GroupElement

export AbstractPermutation, degree, firstmoved, cycle_decomposition

"""
    AbstractPermutation
Abstract type representing permutations of set `1:n`.

Subtypes `Perm <: AbstractPermutation` must implement the following functions:
* `(σ::Perm)(i::Integer)` - the image of `i` under `σ`,
* `degree(σ::Perm)` - the minimal `n` such that `σ(k) == k` for all `k > n`,
* `Perm(images::AbstractVector{<:Integer}[, check::Bool=true])` - construct a
`Perm` from a vector of images. Optionally the second argument `check` may be
set to `false` when the caller knows that `images` constitute a honest
permutation.
"""
abstract type AbstractPermutation <: GroupElement end

"""
    degree(σ::AbstractPermutation)
Return a minimal number `n` such that `σ(k) == k` for all `k > n`.

Such number `n` can be understood as a _degree_ of a permutation, since we can regard
`σ` as an element of `Sym(n)` (and not of `Sym(n-1)`).
By convention `degree` of the trivial permutation must return `1`.
"""
function degree end

"""
    firstmoved(σ::AbstractPermutation)
Return the minimal `i` such that `σ(i) ≠ i`, or `nothing` if `σ == one(σ)`.
"""
firstmoved(σ::AbstractPermutation) = nextmoved(σ, 1)

"""
    nextmoved(σ::AbstractPermutation, d::Integer)
Return the first `i ≥ d` such that `σ(i) ≠ i`, or `nothing` if no such `i` exists.
"""
function nextmoved(σ::AbstractPermutation, d::Integer)
    for i in d:degree(σ)
        if σ(i) ≠ i
            return i
        end
    end
    return nothing # when σ fixes all points ≥ d
end

Base.one(σ::P) where P<:AbstractPermutation = P(Int[], false)
Base.isone(σ::AbstractPermutation) = degree(σ) == 1

function Base.inv(σ::P) where P<:AbstractPermutation
    img = collect(1:degree(σ))
    for i in 1:degree(σ)
        img[σ(i)] = i
    end
    return P(img, false)
end

function Base.:(*)(σ::P, τ::AbstractPermutation) where P<:AbstractPermutation
    deg = max(degree(σ), degree(τ))
    img = collect(1:deg)
    for i in 1:deg
        img[i] = τ(σ(i))
    end
    return P(img, false)
end

function power_by_composition(σ::AbstractPermutation, n::Integer)
    @assert n ≥ 0
    img = collect(1:degree(σ))
    for i in 1:degree(σ)
        k = i
        for j in 1:n
            k = σ(k)
        end
        img[i] = k # = σ(σ(...σ(i)...))
    end
    return typeof(σ)(img, false)
end

function Base.:(^)(σ::AbstractPermutation, n::Integer)
    n < 0 && return inv(σ)^n
    if 1 ≤ n ≤ 10
        return power_by_composition(σ, n)
    else
        # TODO: use cycle decomposition to speedup powering?
        return Base.power_by_squaring(σ, n)
    end
end

function Base.:(==)(σ::AbstractPermutation, τ::AbstractPermutation)
    degree(σ) ≠ degree(τ) && return false
    for i in 1:degree(σ)
        if σ(i) != τ(i)
            return false
        end
    end
    return true
end

Base.hash(σ::AbstractPermutation, h::UInt) =
    foldl((h,i) -> hash(σ(i), h), 1:degree(σ), init=hash(typeof(σ), h))
#=
# has is roughly equivalent to this
function Base.hash(σ::AbstractPermutation, h::UInt)
    h = hash(typeof(σ), h)
    for i in 1:degree(σ)
        h = hash(σ(i), h)
    end
    return h
end
=#

function Base.show(io::IO, σ::AbstractPermutation)
    if isone(σ)
        print(io, "()")
    else
        for cycle in cycle_decomposition(σ)
            if length(cycle) == 1
                continue
            else
                print(io, "(")
                join(io, cycle, ",")
                print(io, ")")
            end
        end
    end
end

function cycle_decomposition(σ::AbstractPermutation)
    isone(σ) && return Vector{Vector{Int}}()
    visited = falses(degree(σ))
    cycles = Vector{Vector{Int}}()
    # each cycle will be a Vector{Int} and we have a whole bunch of them
    for i in 1:degree(σ)
        if visited[i]
            # if we have already seen this point there is no point in computing
            # the same orbit twice
            continue # i.e. skip the rest of the body and continue with the next i
        end
        Δ = _orbit_plain(i, σ, ^)
        visited[Δ] .= true # modify the `visited` along the whole orbit
        push!(cycles, Δ) # add obtained orbit to cycles
    end
    return cycles
end

Base.:^(i::Integer, σ::AbstractPermutation) = σ(i)

function _orbit_plain(pt::Integer, σ::AbstractPermutation, action=^)
    Δ = [pt]
    Δ_set = Set(Δ)
    for δ in Δ
        γ = action(δ, σ)
        if γ ∉ Δ_set
            push!(Δ, γ)
            push!(Δ_set, γ)
        end
    end
    return Δ
end

end # of module AbstractPermutations
