@testset "Image basis" begin

    G = PermutationGroup([perm"(1,2)", perm"(1,2,3)"])

    sc = StabilizerChain(gens(G))
    b = CGT_UniHeidelberg_2022.basis(sc)

    @test CGT_UniHeidelberg_2022.perm_by_images(b, sc) == one(G)
    img = [2,3,1]
    g = CGT_UniHeidelberg_2022.perm_by_images(img, sc)

    @test all(b[i]^g == img[i] for i in 1:length(img))

    w = rand(G, 1000*order(Int, G))
    uniform = length(w)/order(Int, G)
    # we get approximately uniform distribution
    @test abs(count(==(rand(G)), w) - uniform) < 1e-2
end

@testset "Iteration and backtrack" begin

    G = PermutationGroup([perm"(1,2)", perm"(1,2,3)"])
    @test order(G) == 6
    @test Set(collect(G)) == Set([
        perm"(1)",
        perm"(1,2)",
        perm"(1,2,3)",
        perm"(2,3)",
        perm"(1,3,2)",
        perm"(1,3)",
    ])

    using Random
    G, g, h = let n = 12
        g1 = Permutation(Random.randperm(n))
        g2 = Permutation(Random.randperm(n))
        G = PermutationGroup([g1, g2]) # with high probability Sym(12)

        secret = prod(rand(gens(G), 20)) # a random element from G

        g = g1
        h = inv(secret)*g*secret

        G, g, h
    end

    # Feel free to write these tests as you please
    # I'm not forcing a particular syntax for "oracles" and backtrack search

    # Exercise:
    # Find an element in G which sends g to h
end

@testset "Puzzles" begin

    include("puzzles.jl")

    # Exercise: write a function
    # CGT_UniHeidelberg_2022.decompose(
    #     g::AbstractPermutation,
    #     G::AbstractPermGroup,
    # )

    # which decomposes the given element into a product of generators of G
    # to satisfy the tests below:
    let G = PermutationGroup(square_medium)
        S = union(gens(G), inv.(gens(G))) # symmetric generating set
        x = prod(rand(S, 50)) # random element from G

        elts = CGT_UniHeidelberg_2022.decompose(x, G)
        @test elts isa AbstractVector{<:CGT.AbstractPermutation}
        @test all(g ∈ S for g in elts)
        @test isone(x*inv(prod(elts)))
    end

    let G = PermutationGroup(square_hard)
        S = union(gens(G), inv.(gens(G))) # symmetric generating set
        x = prod(rand(S, 50)) # random element from G

        elts = CGT_UniHeidelberg_2022.decompose(x, G)
        @test elts isa AbstractVector{<:CGT.AbstractPermutation}
        @test all(g ∈ S for g in elts)
        @test isone(x*inv(prod(elts)))
    end

    # Exercise: use the backtrack search to bring a randomly scrambled cube to
    # its solving position, i.e. every element of the top face is on the top
    # face (possibly in arbitrary order), etc.
    # This amounts to setting an oracle which accepts every permutation which maps
    # {1..8 }->{1..8}
    # {9..16}->{9..16}, etc.


    # let G = PermutationGroup(cube4)
    #     S = union(gens(G), inv.(gens(G)))
    #     x = prod(rand(S, 50))

    #     g = # backtrack here on stabilizer chain of G

    #     elts = CGT_UniHeidelberg_2022.decompose(g, G)
    #     @test elts isa AbstractVector{<:AbstractPermutation}
    #     @test all(g ∈ S for g in elts)
    #     t = x*inv(prod(elts))
    #     for face in (1:8, 9:16, 17:24, 25:32, 33:40, 41:48)
    #         F = BitSet(face)
    #         @test all(i^t ∈ F for i in face)
    #     end
    # end
end
