# included from CGT_UniHeidelberg_2022.jl

# Some possible methods:
# * Is a word a generator or a relation (this does not necessarily have to be part of the struct, to save space)
# * Is a word the identity (word problem) or reduced
# * More generally, are two words in the same equivalence class (w ∼ z)
# * Group multiplication / concatenation
# * (Optional) reference to Alphabet, e.g. for word comparison
abstract type AbstractWord{T} <: AbstractArray{T, 1} end

export Word

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

# AbstractArray interface
Base.size(w::Word) = size(w.letter_indices)
Base.getindex(w::Word, i::Int) = getindex(w.letter_indices[i])

function Base.setindex!(w::Word, letter, i::Int)
    w.letter_indices[i] = letter
end

# Return "mutable array"
function Base.similar(w::Word{T}) where T
    return Word{T}(Vector{T}(undef, length(w)))
end

function Base.similar(::Word{T}, dim::Int) where T
    return Word{T}(Vector{T}(undef, dim))
end

Base.resize!(w::Word, n::Int) = resize!(w.letter_indices, n)

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

# The identity is represented by the empty word
Base.one(::Word{T}) where T = Word{T}(T[])
Base.isone(w::AbstractWord) = isempty(w)

# Taken from 05_Alphabets_and_words.jl
function inv!(out::AbstractWord, A::Alphabet, w::AbstractWord)
    @assert length(out) == length(w)
    for (idx, letter) in enumerate(Iterators.reverse(w))
	out[idx] = inv(A, letter)
    end
    return out
end
Base.inv(A::Alphabet, w::AbstractWord) = inv!(similar(w), A, w)

# Serialization
function Base.show(io::IO, ::MIME"text/plain", w::AbstractWord{T}) where T
    l = length(w)
    for (i, letter) in enumerate(w)
	    print(io, letter)
	    if i < l
	        print(io, '·')
	    end
    end
end

string_repr(A::Alphabet, w::AbstractWord) = join((A[idx] for idx in w), '·')

# When rewriting a word to reduced form, we can use a queue approach
# where letters followed by their inverse are removed. This is done by
# exposing a push! and pop!-style functions.
Base.push!(w::Word, letter_idx) = push!(w.letter_indices, letter_idx)
Base.pushfirst!(w::Word, letter_idx) = pushfirst!(w.letter_indices, letter_idx)
Base.pop!(w::Word) = pop!(w.letter_indices)
Base.popfirst!(w::Word) = popfirst!(w.letter_indices)
Base.append!(w::Word, letter_idx) = append!(w.letter_indices, letter_idx)
Base.prepend!(w::Word, letter_idx) = prepend!(w.letter_indices, letter_idx)

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
