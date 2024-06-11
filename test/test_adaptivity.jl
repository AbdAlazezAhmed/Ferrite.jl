@testset "Kopp Quadrilateral Adaptivity" begin
    @testset "Kopp Neighborhood" begin
        grid = generate_grid(Ferrite.KoppCell, (2,2))
        # No refinement
        @test Ferrite.get_neighbors(grid, 2=>0) == (-1=>0, -1=>0, 4=>0, 1=>0)
        # Refine one level
        Ferrite.refine!(grid, 1=>0)
        @test Ferrite.get_neighbors(grid, 2=>0) == (-1=>0, -1=>0, 4=>0, 1=>0)
        @test Ferrite.get_neighbors(grid, 1=>3) == (1=>2, 2=>0, 3=>0, 1=>4)
        Ferrite.refine!(grid, 2=>0)
        @test Ferrite.get_neighbors(grid, 2=>1) == (-1=>0, 2=>2, 2=>4, 1=>2)
        @test Ferrite.get_neighbors(grid, 1=>3) == (1=>2, 2=>4, 3=>0, 1=>4)
        # Refine two levels
        Ferrite.refine!(grid, 1=>3)
        @test Ferrite.get_neighbors(grid, 1=>6) == (1=>2, 2=>4, 1=>7, 1=>5)
        @test Ferrite.get_neighbors(grid, 1=>7) == (1=>6, 2=>4, 3=>0, 1=>8)
        @test Ferrite.get_neighbors(grid, 2=>4) == (2=>1, 2=>3, 4=>0, 1=>3)
        Ferrite.refine!(grid, 2=>4)
        @test Ferrite.get_neighbors(grid, 2=>5) == (2=>1, 2=>6, 2=>8, 1=>6)
        @test Ferrite.get_neighbors(grid, 2=>8) == (2=>5, 2=>7, 4=>0, 1=>7)
        @test Ferrite.get_neighbors(grid, 1=>6) == (1=>2, 2=>5, 1=>7, 1=>5)
        @test Ferrite.get_neighbors(grid, 1=>7) == (1=>6, 2=>8, 3=>0, 1=>8)
    end

    @testset "Kopp Dofs" begin
        grid = generate_grid(Ferrite.KoppCell, (2,2))
        ip = Lagrange{RefQuadrilateral, 1}()
        qr = QuadratureRule{RefQuadrilateral}(2)
        cellvalues = CellValues(qr, ip); 
        @testset "Single SubDofHandler" begin
            dh = DofHandler(grid)
            add!(dh, :u, ip)
            close!(dh);
        end
    end

end