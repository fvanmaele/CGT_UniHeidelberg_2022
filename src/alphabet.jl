# included from CGT_UniHeidelberg_2022.jl

export Alphabet, hasinverse, setinverse!

""" Allow arbitrary objects (except Integers) as letters
Types:
 * `T`: object representing a letter in the alphabet
Indices are assumed to be of type Int64, the same as with Vector{T}.
"""
struct Alphabet{T}
    gens::Vector{T}
    gens_to_int::Dict{T, Int}  # maps objects T to integers
    inv_state::Dict{T,T}       # maps objects T to formal inverses

    function Alphabet{T}(S::AbstractVector) where T
        # To distinguish between "indexing integers" and "indexing letters" we
        # do not allow integers as element types.
        @assert !(T <: Integer)

        # The set may either contain redundant, or self-inverse generators. In
        # this case, integer representations are not unique (with the last
        # occurence defining the index).  To avoid this case, verify that the
        # set of generators has no redundant elements.
        @assert (length(unique(S)) == length(S))
        A = new{T}(S, Dict{T,typeof(1)}(), Dict{T,T}())

        # Map integers sequentially (starting at 1)
        for i in 1:length(S)
            push!(A.gens_to_int, S[i]=>i)
        end
        return A
    end
end

# Loop methods
Base.iterate(A::Alphabet) = iterate(A.gens)
Base.iterate(A::Alphabet, state) = iterate(A.gens, state)
Base.length(A::Alphabet) = length(A.gens)

# O(1) element check
Base.in(letter, A::Alphabet) = haskey(A.gens_to_int, letter)

# Maps an index to a letter
function Base.getindex(A::Alphabet, index::Integer)
    @assert index >= 1 "alphabet indices start at 1"
    @assert index <= length(A.gens) "alphabet is shorter than requested $index"
    return A.gens[index]
end

# Maps a letter to an index
function Base.getindex(A::Alphabet, letter)
    @assert letter ∈ A "letter is not in alphabet"
    return A.gens_to_int[letter]
end

function setinverse!(A::Alphabet{T}, x::T, X::T) where T
    @assert x ∈ A "letter is not in alphabet"
    @assert X ∈ A "letter is not in alphabet"
    A.inv_state[x] = X
    A.inv_state[X] = x
end

# Q: The exercise mentions when asking for an inverse, and the element is not
# invertible, that an error should be returned. What would be the benefit in
# this case over return a boolean? (An error is raised in the `inv` function)
function hasinverse(A::Alphabet, letter)
    @assert letter ∈ A "letter is not in alphabet"
    return haskey(A.inv_state, letter)
end

function hasinverse(A::Alphabet, index::Integer)
    return haskey(A.inv_state, A[index])
end

function Base.inv(A::Alphabet, letter)
    if hasinverse(A, letter)
        return A.inv_state[letter] # get inverted letter
    else
        throw(DomainError("letter has no inverse"))
        #return nothing
    end
end

function Base.inv(A::Alphabet, index::Integer)
    letter = A[index]
    if hasinverse(A, letter)
        return A[A.inv_state[letter]] # get index of inverted letter
    else
        throw(DomainError("letter has no inverse"))
        #return nothing
    end
end

# Note: only prints indices, not elements itself
function Base.show(io::IO, A::Alphabet{T}) where T
    println(io, "Alphabet of $T with $(length(A)) letters:")
    for letter in A
        print(io, A[letter], "\t", letter)
        if hasinverse(A, letter)
            println(io, " with inverse ", A[inv(A, A[letter])])
        else
            println(io, "")
        end
    end
end
