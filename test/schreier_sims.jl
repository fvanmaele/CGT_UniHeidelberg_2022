@testset "Schreier-Sims algorithm" begin
    order(::Type{I}, sc::CGT.StabilizerChain) where I<:Integer =
        prod(I∘length, CGT.transversals(sc))
    order(sc::CGT.StabilizerChain) = order(BigInt, sc)

    # Make sure the order of generators does not matter
    S1 = [perm"(1,2)", perm"(1,2,3)"]
    sc1 = CGT.schreier_sims(S1)
    S2 = [perm"(1,2,3)", perm"(1,2)"]
    sc2 = CGT.schreier_sims(S2)
    @test order(sc1) == order(sc2)
    @test order(sc1) == 6

    for group_order in 2:30
        for S in SmallPermGroups[group_order]
            sc = CGT.schreier_sims(S) # defaults to Transversal
            @test order(Int, sc) == group_order
        end
    end

    for group_order in 2:30
        for S in SmallPermGroups[group_order]
            STree_t = CGT.SchreierTree{Int, eltype(S), typeof(^)}
            sc_tree = CGT.schreier_sims(STree_t, S)
            @test order(Int, sc_tree) == group_order
        end
    end
end
