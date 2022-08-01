using Test
using CGT_UniHeidelberg_2022
const CGT = CGT_UniHeidelberg_2022

import CGT_UniHeidelberg_2022:
    Permutation,
    CyclePermutation,
    CyclePermutation2,
    Orbit,
    Transversal,
    Alphabet

include("small_perm_groups.jl")

@testset verbose=true "CGT" begin
    include("permutations.jl")
    include("transversals.jl")
    include("schreier_sims.jl")
    include("alphabet.jl")
    include("word.jl")
end
