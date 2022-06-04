"""
    AbstractOrbit{T}
Represents an orbit of elements of type `T` under an action of a group from the right.

Subtypes of `AbtractOrbits` must define functions from iterator interface i.e.
* `Base.iterate` and `Base.length`;

Additionally the following need to be implemented
* `Base.in` - determines if a point is in the orbit
* `Base.first` - return the distinguished element of the orbit
* `action` - return the action used to create the orbit
* `Base.push!` - add a point to the orbit

`Base.push!` should be a no-op (trivial) if pushed point is already in the orbit.
"""
abstract type AbstractOrbit{T} end

Base.eltype(::Type{<:AbstractOrbit{T}}) where T = T

struct Orbit{T,A} <: AbstractOrbit{T}
    points::Vector{T}
    set::Set{T}
    action::A

    function Orbit{T, A}(pts::AbstractVector, S::AbstractVector{<:GroupElement}, action=^) where {T,A}
        @assert !isempty(S)
        Δ = new{T, A}(pts, Set(pts), action)
        for δ in Δ
            for s in S
                push!(Δ, action(δ, s))
            end
        end
        return Δ
    end
end

Orbit(pts::AbstractVector, S::AbstractVector{<:GroupElement}, action=^) =
    Orbit{eltype(pts), typeof(action)}(pts, S, action)
Orbit(pt, S::AbstractVector{<:GroupElement}, action=^) = Orbit([pt], S, action)

Base.length(o::Orbit) = length(o.points)
Base.iterate(o::Orbit) = iterate(o.points)
Base.iterate(o::Orbit, state) = iterate(o.points, state)
Base.in(x, o::Orbit) = x in o.set
Base.first(o::Orbit) = first(o.points)
action(o::Orbit) = o.action

function Base.push!(o::Orbit, x)
    x in o && return o

    push!(o.points, x)
    push!(o.set, x)
    return o
end

"""
    AbstractTransversal{T, GEl<:GroupElement}
Transversals are data structures representing orbits which additionally carry
a group representative for its elements.

Additionally to the `AbstractOrbit` methods `AbstractTransversal`s must
implement a method to retrive the associated group element through
 * `Base.getindex(T::AbstractTransversal{T}, y)` which must satisfy
 `T[y] = g` whenever `action(T)(first(T), g) = y` and throw a `KeyError` if
 `y` is not in orbit.
 * `Base.push!` with `y=>g` pair which stores representative `g` together with
 point `y`. Note: this operation should replace the existing representative if
 `y` is already present in the transversal.
"""
abstract type AbstractTransversal{T,GEl<:GroupElement} <: AbstractOrbit{T} end

struct Transversal{T,GEl,Ac} <: AbstractTransversal{T,GEl}
    points::Vector{T}
    representatives::Dict{T,GEl}
    action::Ac

    function Transversal{T,GEl,Ac}(
        point,
        S::AbstractVector{<:GroupElement},
        action = ^,
    ) where {T,GEl,Ac}
        @assert !isempty(S) "The generating set must be non-empty, got $S"
        trans = new{T,GEl,Ac}(T[point], Dict(point => one(first(S))), action)
        for pt in trans
            for s in S
                y = action(pt, s)
                y ∈ trans && continue
                push!(trans, y=>trans[pt]*s)
            end
        end
        return trans
    end
end

Transversal(point::T, S::AbstractVector{GEl}, action = ^) where {T, GEl<:GroupElement}=
    Transversal{T,GEl,typeof(action)}(point, S, action)

Base.iterate(t::Transversal) = iterate(t.points)
Base.iterate(t::Transversal, state) = iterate(t.points, state)
Base.length(t::Transversal) = length(t.points)
Base.in(x, t::Transversal) = haskey(t.representatives, x)
Base.first(t::Transversal) = first(t.points)
action(t::Transversal) = t.action

function Base.getindex(t::Transversal, pt)
    pt ∉ t && throw(KeyError(pt))
    return t.representatives[pt]
end

function Base.push!(t::Transversal, y_g::Pair{T, <:GroupElement}) where T
    y, g = y_g
    if !(y in t)
        push!(t.points, y)
    end
    t.representatives[y] = g
    return t
end
Base.setindex!(t::Transversal, g::GroupElement, pt) = push!(t, pt=>g)

# Assignment 3
# Implement struct SchreierTree{T, GEl, ...} <: AbstractTransversal{T, GEl} ... end
struct SchreierTree{T,GEl,Ac} <: AbstractTransversal{T,GEl}
    points::Vector{T}
    factors::Dict{T,GEl} # factors (generators) for representatives in transversal
    action::Ac

    function SchreierTree{T,GEl,Ac}(
        point,
        S::AbstractVector{<:GroupElement},
        action = ^,
    ) where {T,GEl,Ac}
        @assert !isempty(S) "The generating set must be non-empty, got $S"
        Sch = new{T,GEl,Ac}(T[point], Dict(point => one(first(S))), action)
        for pt in Sch
            for s in S
                y = action(pt, s)
                # Q: The AbstractTransversal interface dicates that "[push!] should replace the 
                # existing representative if `y` is already present in the transversal.". This is
                # done accordingly in the implementation of push!, which includes an `in` check.
                # However, this `in` check is repeated here, with a `continue` clause, preventing
                # the new representative to be pushed. Uncommenting the line, however, results in
                # a non-trivial factor for the distinguished element, and getindex does
                y ∈ Sch && continue
                push!(Sch, y=>s) # Push (reference to) group generator
            end
        end
        return Sch
    end
end

SchreierTree(point::T, S::AbstractVector{GEl}, action = ^) where {T, GEl<:GroupElement}=
    SchreierTree{T,GEl,typeof(action)}(point, S, action)

Base.iterate(t::SchreierTree) = iterate(t.points)
Base.iterate(t::SchreierTree, state) = iterate(t.points, state)
Base.length(t::SchreierTree) = length(t.points)
Base.in(x, t::SchreierTree) = haskey(t.factors, x)
Base.first(t::SchreierTree) = first(t.points)
action(t::SchreierTree) = t.action

# Begin at last factor, and work our way up (code adapted from Assignment 2)
function Base.getindex(t::SchreierTree, pt)
    pt ∉ t && throw(KeyError(pt))
    current = pt
    g = one(t.factors[current])

    while (gen = t.factors[current]) != one(t.factors[current])
        current = t.action(current, inv(gen))
        g = gen * g
    end
    return g
end

# Note that we push generators here, instead of group elements.
function Base.push!(t::SchreierTree, y_s::Pair{T, <:GroupElement}) where T
    y, s = y_s
    if !(y in t)
        push!(t.points, y)
    end
    t.factors[y] = s
    return t
end
Base.setindex!(t::SchreierTree, g::GroupElement, pt) = push!(t, pt=>g)
