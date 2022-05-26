include("../src/collections.jl")
using Test
#
#using Main.sample
#
#Setting up some collections
col312_1 = [[1, 2, 4], [4, 5, 7], [7, 8, 10], [1, 10, 11], [1, 2, 5], [4, 5, 8], [7, 8, 11], [2, 10, 11], [1, 2, 11], [2, 4, 5], [5, 7, 8], [8, 10, 11], [2, 5, 8], [5, 8, 11], [2, 8, 11], [2, 5, 11], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7], [6, 7, 8], [7, 8, 9], [8, 9, 10], [9, 10, 11], [10, 11, 12], [1, 11, 12], [1, 2, 12]]
c1 = LabelCollection(3, 12, col312_1)
col312_2 = [[1, 2, 4], [4, 5, 7], [7, 8, 10], [1, 10, 11], [1, 2, 10], [1, 4, 5], [4, 7, 8], [7, 10, 11], [1, 2, 11], [2, 4, 5], [5, 7, 8], [8, 10, 11], [1, 4, 7], [4, 7, 10], [1, 7, 10], [1, 4, 10], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7], [6, 7, 8], [7, 8, 9], [8, 9, 10], [9, 10, 11], [10, 11, 12], [1, 11, 12], [1, 2, 12]]
c2 = LabelCollection(3, 12, col312_2)

NC = [[4,5,8], [3, 5, 7], [7, 8, 10], [1, 10, 11], [1, 2, 5], [4, 5, 8], [7, 8, 11], [2, 10, 11], [1, 2, 11], [2, 4, 5], [5, 7, 8], [8, 10, 11], [2, 5, 8], [5, 8, 11], [2, 8, 11], [2, 5, 11], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7], [6, 7, 8], [7, 8, 9], [8, 9, 10], [9, 10, 11], [10, 11, 12], [1, 11, 12], [1, 2, 12]]
NCL = LabelCollection(3, 12, NC)

@testset "Collections" begin
    @testset "NonCrossing" begin
        # Testing a non crossing collection
        @test !isNonCrossingCollection(NCL);

        #Noncrossing
        @test isNonCrossingCollection(c1);
        @test isNonCrossingCollection(c2);
    end

    @testset "Maximal Collection" begin
        @test isMaximalNonCrossingCollection(c1)
        @test isMaximalNonCrossingCollection(c2)
        @test !isMaximalNonCrossingCollection(NCL)
    end
end
