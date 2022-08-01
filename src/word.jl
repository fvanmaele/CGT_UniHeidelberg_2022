# included from CGT_UniHeidelberg_2022.jl

export AbstractWord, Word

# AbstractWord API imported from notebooks/06_Rewriting.jl
"""
    AbstractWord{T} <: AbstractVector{T}
Abstract type representing words over an Alphabet.

`AbstractWord` is just a string of integers and as such gains its meaning in the
contex of an Alphabet (when integers are understood as pointers to letters).
The subtypes of `AbstractWord{T}` need to implement the following methods which
constitute `AbstractWord` interface:
 * a constructor from `AbstractVector{T}`
 * linear indexing (1-based) consistent with iteration returning pointers to letters of an alphabet (`getindex`, `setindex!`, `size`),
 * `Base.push!`/`Base.pushfirst!`: append a single value at the end/beginning,
 * `Base.pop!`/`Base.popfirst!`: pop a single value from the end/beginning,
 * `Base.append!`/`Base.prepend!`: append a another word at the end/beginning,
 * `Base.resize!`: drop/extend a word at the end to the requested length
 * `Base.similar`: an uninitialized word of a similar type/storage.

Note that `length` represents free word (how it is written in an alphabet)
and not its the shortest form (e.g. the normal form).

!!! note
    It is assumed that `eachindex(w::AbstractWord)` returns `Base.OneTo(length(w))`

The following are implemented for `AbstractWords` but can be overloaded for
performance reasons:

* `Base.==`: the equality (as words),
* `Base.hash`: simple uniqueness hashing function
* `Base.:*`: word concatenation (monoid binary operation),
"""
abstract type AbstractWord{T} <: AbstractVector{T} end

function Base.hash(w::AbstractWord, h::UInt)
    return foldl((h, x) -> hash(x, h), w, init = hash(AbstractWord, h))
end

@inline function Base.:(==)(w::AbstractWord, v::AbstractWord)
    length(w) == length(v) || return false
    return all(w[i] == v[i] for i in eachindex(w))
end

# resize! + copyto!
function store!(w::AbstractWord, v::AbstractWord)
    resize!(w, length(v))
    copyto!(w, v)
    return w
end

# The identity is represented by the empty word
Base.one(::Type{W}) where {T,W<:AbstractWord{T}} = W(T[])
Base.one(::W) where {W<:AbstractWord} = one(W)
Base.isone(w::AbstractWord) = isempty(w)

function Base.getindex(w::W, u::AbstractRange) where {W<:AbstractWord}
    return W([w[i] for i in u])
end

function Base.:^(w::AbstractWord, n::Integer)
    return n >= 0 ? repeat(w, n) :
           throw(
        DomainError(
            n,
            "To rise a Word to negative power you need to provide its inverse.",
        ),
    )
end
	
function Base.literal_pow(::typeof(^), w::AbstractWord, ::Val{p}) where {p}
    return p >= 0 ? repeat(w, n) :
           throw(
        DomainError(
            p,
            "To rise a Word to negative power you need to provide its inverse.",
        ),
    )
end
	
function Base.show(io::IO, ::MIME"text/plain", w::AbstractWord)
    print(io, typeof(w), ": ")
    return show(io, w)
end

function Base.show(io::IO, w::AbstractWord{T}) where {T}
    if isone(w)
        print(io, "(id)")
    else
        join(io, w, "·")
    end
end


""" Word structure with flexible storage type (default to UInt8)
"""
struct Word{T} <: AbstractWord{T}
    letter_indices::Vector{T}

    function Word{T}(idx::AbstractVector) where T
        w = new(idx)
        return w
    end
end

Word(idx::AbstractVector) = Word{UInt8}(idx)

# AbstractVector interface
Base.size(w::Word) = size(w.letter_indices)
Base.getindex(w::Word, i::Integer) = w.letter_indices[i]
Base.setindex!(w::Word, v, i::Integer) = w.letter_indices[i] = v

Base.resize!(w::Word, n::Integer) = resize!(w.letter_indices, n)
Base.similar(w::Word{T}) where T = Word{T}(similar(w.letter_indices))
Base.similar(w::Word{T}, n::Integer) where T = Word{T}(similar(w.letter_indices, n))

# * multiplication
# When concatenating words, we can either proceed naively and
# concatenate both words, or we can try to make reduction in the
# process, i.e. by the rules:
#  - ysŝx ∼ yx
#  - yŝsx ∼ yx
# In the former case, we do not need an Alphabet, but the length of
# `out` is not known a-priori.
function mul!(out::AbstractWord, w::AbstractWord, z::AbstractWord)
    resize!(out, length(w)+length(z))
    copyto!(out, 1, w, 1, length(w))
    copyto!(out, length(w)+1, z, 1, length(z))
    return out
end

function Base.:*(w::AbstractWord, z::AbstractWord)
    out = similar(w, length(w) + length(z))
    mul!(out, w, z)
    return out
end

# Taken from 05_Alphabets_and_words.jl
function inv!(out::AbstractWord, A::Alphabet, w::AbstractWord)
    @assert length(out) == length(w)
    for (idx, letter) in enumerate(Iterators.reverse(w))
	    out[idx] = inv(A, letter)
    end
    return out
end
Base.inv(A::Alphabet, w::AbstractWord) = inv!(similar(w), A, w)

# When rewriting a word to reduced form, we can use a queue approach
# where letters followed by their inverse are removed. This is done by
# exposing a push! and pop!-style functions.
Base.pop!(w::Word) = pop!(w.letter_indices)
Base.popfirst!(w::Word) = popfirst!(w.letter_indices)

Base.push!(w::Word, n::Integer) = push!(w.letter_indices, n)
Base.pushfirst!(w::Word, n::Integer) = pushfirst!(w.letter_indices, n)

Base.append!(w::Word, v::AbstractWord) = append!(w.letter_indices, v)
Base.prepend!(w::Word, v::AbstractWord) = prepend!(w.letter_indices, v)

# Rewriting rules:
# * xsŝy → xy  [1]
# * xŝsy → xy  [2]
function Base.reduce(A::Alphabet, w::Word)
    queue = similar(w, 0)
    cnt = 1

    # go over each element of the input word
    while cnt <= length(w)
        idx1 = w[cnt]
        cnt += 1

        if length(queue) > 0
            idx2 = queue[length(queue)] # previous element
            l1 = A[idx1]  # returns Alphabet[index] → letter{T}
            l2 = A[idx2]

            if (l1 == inv(A, l2)) || (l2 == inv(A, l1))
                # remove previous letter from input queue
                pop!(queue)
            else
                # otherwise, push the new letter for processing
                push!(queue, idx1)
            end
        else
            push!(queue, idx1)
        end
    end
    return queue
end

# BufferWord (deque)
