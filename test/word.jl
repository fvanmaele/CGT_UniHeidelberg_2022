@testset "Word API" begin
    # Construct alphabet
    σ = Permutation([2,1,3])
    τ = Permutation([1,3,2])
    A = Alphabet{Permutation}([σ, τ])
    setinverse!(A, σ, σ)
    setinverse!(A, τ, τ)

    # Word indices
    idx1 = [1, 1, 2, 1, 2, 2]
    idx2 = [2, 1]
    idx3 = [2]
    idx4 = [1, 1, 2, 1, 1, 2]

    W1 = Word(idx1);
    W2 = Word(idx2);
    W3 = Word(idx3);
    W4 = Word(idx4);

    # Check default type is UInt8
    @test eltype(W1) == UInt8;

    @testset "Word operations" begin
        # Multiplication
        W11 = W1 * W1;
        @test length(W11) == length(W1)*2
        @test collect(W11) == [collect(W1); collect(W1)]
        @test W11 == repeat(W1, 2) # this guy is a `Word`!
        @test isempty(one(W1))

        # Inverses
        # Note: elements in A are self-inverse
        @test inv(A, W1) == Word([2, 2, 1, 2, 1, 1])
        @test inv(A, W2) == Word([1, 2])
        @test inv(A, W3) == Word([2])
        @test inv(A, W4) == Word([2, 1, 1, 2, 1, 1])

        # Reduction
        @test reduce(A, W1) == Word([2, 1])
        @test reduce(A, W2) == Word([2, 1])
        @test reduce(A, W3) == Word([2])
        @test reduce(A, W4) == one(W4)

        @test reduce(A, W1 * inv(A, W1)) == one(W1)
        @test reduce(A, W2 * inv(A, W2)) == one(W2)
        @test reduce(A, W3 * inv(A, W3)) == one(W3)
        @test reduce(A, W4 * inv(A, W4)) == one(W4)
    end
end
