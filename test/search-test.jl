@testset "multiplication_search k = GF(3^2)..." begin

    k, x = Nemo.FiniteField(3, 2, "x")
    A = multiplication_search(k, (0, 0), 3)
    @test length(A) == 1
    @test A[1] == ((1, 2), (2*x+1, 2), (x, 2))

end

@testset "multiplication_search k = GF(3^3)..." begin

    k, x = Nemo.FiniteField(3, 3, "x")
    A = multiplication_search(k, (0, 0, 0), 6)
    @test length(A) == 1
    @test A[1] == ((x^2+1, 1), (x^2+2*x+1, 2), (2*x^2+x+1, 1), (x^2+x, 2), (2*x^2+x, 2), (x^2, 1))

    A = multiplication_search(k, (1, 0, 0), 6)
    @test length(A) == 2
    @test A[1] == ((x^2+1, 1), (x^2+2*x+1, 2), (2*x^2+x+1, 1), (x^2+x, 2), (2*x^2+x, 2), (x^2, 1))
    @test A[2] == ((x^2+1, 2), (x^2+x+1, 1), (2*x^2+1, 2), (2*x^2+2*x+1, 2), (2*x^2+x, 1), (x^2, 2))

    A = multiplication_search(k, (2, 1, 0), 6)
    @test length(A) == 4
    @test A == [((x+1, 1), (2*x+1, 2), (2*x^2+1, 1), (x, 1), (x^2+x, 2), (2*x^2+x, 1)),
 ((x^2+1, 1), (x^2+2*x+1, 2), (2*x^2+x+1, 1), (x^2+x, 2), (2*x^2+x, 2), (x^2, 1)),
 ((x^2+1, 2), (x^2+x+1, 1), (2*x^2+1, 2), (2*x^2+2*x+1, 2), (2*x^2+x, 1), (x^2, 2)),
 ((x^2+x+1, 2), (x^2+2*x+1, 1), (2*x^2+1, 1), (2*x^2+x+1, 2), (2*x^2+2*x+1, 1), (x^2+x, 1))]

    A = multiplication_search(k, (2, 1, 0), 5)
    @test length(A) == 0
end
