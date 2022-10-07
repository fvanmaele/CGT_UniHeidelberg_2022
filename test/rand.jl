@testset "Random element generation" begin
    # S = SmallPermGroups[10][2]  # order 10
    # sc = CGT.schreier_sims(S)   # depth 2

    for group_order in 2:30
        for S in SmallPermGroups[group_order]
            sc = CGT.schreier_sims(S)

            # check that elements are (approximately) uniformly distributed
            V = []
            n_reps = 50000
            for i = 1:n_reps
                push!(V, rand(sc))
            end

            # get counts
            C = [(i, count(==(i), V)) for i in unique(V)]
            @test length(C) == group_order # all elements were hit
            for i = 1:length(C)
                pct = C[i][2] / n_reps
                @test pct > 0.9 * ((n_reps / group_order) / n_reps)  # 90% confidence
                @test pct < 1.1 * ((n_reps / group_order) / n_reps)
            end
        end
    end
end
