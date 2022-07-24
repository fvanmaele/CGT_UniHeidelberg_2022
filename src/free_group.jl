# included from CGT_UniHeidelberg_2022.jl

abstract type AbstractFPGroup{T, W} end

export FreeGroup, FPGroupElement

struct FreeGroup{T, W} <: AbstractFPGroup{T, W}
    A::Alphabet{T}
    gens::Vector{W} # stores `Word`s
end

rank(G::FreeGroup) = length(G.gens)

# One may argue to use a reduced word as (representative for) a free
# group element. From a practical (performance) perspective, I see the
# following downsides:
#
# * Reduction is not lazy, e.g. the free group multiplication
#   [w1]*[w2]*[w3] would result in 5 reductions (once for each word,
#   once for each multiplication) instead of 2;
#
# * It is unclear how to pre-allocate reduced words, since the
#   length after reduction is a-priori not known.
#
# As such this implementation stores an arbitrary representative (word)
# of a free group element. Reduction can be done on-demand with the
# mulreduce! and reduce! methods.
#
# Q: If a FPGroupElement is `mutable`, how useful are various `out`
# parameters for performance reasons?
mutable struct FPGroupElement{G<:AbstractFPGroup, W<:AbstractWord} <: GroupElement
    parent::G
    word::W
end

# Pre-allocation methods
function Base.similar(g::FPGroupElement)
    return FPGroupElement(g.parent, similar(g.word))
end

function Base.similar(g::FPGroupElement, dim::Int)
    return FPGroupElement(g.parent, similar(g.word, dim))
end

# Accessors
Base.parent(g::FPGroupElement) = g.parent
word(g::FPGroupElement) = g.word
alphabet(G::FreeGroup) = G.A

# Group operations
function mul!(out::FPGroupElement, g::FPGroupElement, h::FPGroupElement)
    @assert parent(g) === parent(h)
    # Note: `word` returns a mutable object bound to `g`, so this is pure sugar
    mul!(word(out), word(g), word(h))
    return out
end

reduce(G::FreeGroup, g::FPGroupElement) = reduce(alphabet(G), word(g))
# XXX: use modifying reduction
function reduce!(g::FPGroupElement)
    g.word = reduce(parent(g), word(g))
    return g.word
end

# Group multiplication with word reduction.
# To allow reusing the temporary output buffer for different words (of
# the same length), remaining entries are set to 0.
# XXX: whether this trick results in better performance is unknown.
function mulreduce!(out::FPGroupElement, g::FPGroupElement, h::FPGroupElement; pad=true)
    @assert parent(g) === parent(h)
    # concatenation
    mul!(word(out), word(g), word(h))

    # reduction (+1 allocation)
    tmp = reduce(g.parent.A, word(out))
    out.word[1:length(tmp)] = tmp

    # optional padding
    if pad == true
        out.word[length(tmp):length(word(out))] .= 0
    end
    return out, length(tmp)
end

# Simple multiplication, with reduction only done lazly/when necessary
function Base.:*(g::FPGroupElement, h::FPGroupElement)
    @assert parent(g) === parent(h)
    out = similar(g, length(word(g)) + length(word(h)))
    out = mul!(out, word(g), word(h))
    return FPGroupElement(out, parent(g))
end

# Group inversion
# No reduction is performed here: if a word is reduced, then its inverse
# is also reduced.
function inv!(out::FPGroupElement, g::FPGroupElement)
    @assert parent(out) === parent(g)
    inv!(word(out), alphabet(parent(g)), word(g))
    return out
end
Base.inv(g::FPGroupElement) = inv!(similar(g), g)

# XXX: if words are often compared for equality, some kind of caching
# could be implemented (e.g. as an additional field in FPGroupElement),
# or when it is known a set of words is reduced, use == on the
# underlying words (g.W) directly.
function Base.:(==)(g::FPGroupElement, h::FPGroupElement)
    @assert parent(g) === parent(h)
    return reduce(alphabet(parent(g)), word(g)) == reduce(alphabet(parent(h)), word(h))
end

# real interesting task:
# how and when to freely reduce the underlying words
