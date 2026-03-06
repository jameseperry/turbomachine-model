@testset "Utility LinearMap" begin
    U = TurboMachineModel.Utility

    map = U.LinearMap([1.0, 2.0, 4.0], [10.0, 20.0, 40.0])
    @test U.linear_evaluate(map, 1.5) ≈ 15.0
    @test U.linear_evaluate(map, 3.0) ≈ 30.0
    @test U.linear_evaluate(map, 0.0) ≈ 10.0
    @test U.linear_evaluate(map, 5.0) ≈ 40.0

    @test U.linear_evaluate([1.0, 2.0, 4.0], [10.0, 20.0, 40.0], 1.5) ≈ 15.0
    @test U.linear_evaluate([1.0, 2.0, 4.0], [10.0, 20.0, 40.0], 3.0) ≈ 30.0

    @test_throws ErrorException U.LinearMap([1.0], [2.0])
    @test_throws ErrorException U.LinearMap([1.0, 1.0], [2.0, 3.0])
    @test_throws ErrorException U.LinearMap([1.0, 2.0], [2.0])
    @test_throws ErrorException U.linear_evaluate([1.0, 2.0], [3.0], 1.5)
end
