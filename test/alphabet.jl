@testset "Alphabet API" begin 
    σ = Permutation([2,1,3])
    τ = Permutation([1,3,2])
    A = Alphabet{Permutation}([σ, τ])
    
    # Cycles are self-inverse. In this case, the Alphabet constructor should
    # throw an exception.
    @test_throws AssertionError Alphabet{Permutation}([σ, inv(σ), τ, inv(τ)])
    # We do not allow alphabets of integers.
    @test_throws AssertionError Alphabet{Int64}([1, 2, 3])

    @test length(A) == 2
    @test A[1] == σ
    @test A[2] == τ
    @test A[σ] == 1
    @test A[τ] == 2

    κ = Permutation([2,1,3,5,4])
    # Invalid element access
    @test_throws AssertionError A[3]
    @test_throws AssertionError A[κ]

    # Inverse is not defined a-priori
    @test !hasinverse(A, σ)
    @test_throws DomainError inv(A, σ)
    @test_throws DomainError inv(A, 1)

    setinverse!(A, σ, σ) # σ is a cycle
    @test hasinverse(A, σ)
    @test inv(A, σ) == σ
    @test inv(A, 1) == A[σ]

    @test !hasinverse(A, τ)
    @test_throws DomainError inv(A, τ)
    @test_throws DomainError inv(A, 2)

    setinverse!(A, τ, τ) # τ is a cycle
    @test hasinverse(A, τ)
    @test inv(A, τ) == inv(τ)
    @test inv(A, A[τ]) == A[inv(τ)]

end
