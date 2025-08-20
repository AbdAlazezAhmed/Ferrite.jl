include("/home/amohamed/.julia/dev/Ferrite/src/Grid/Adaptivity/kopp.jl")
# @testset "Kopp Quadrilateral Refinement" begin
#     ip = DiscontinuousLagrange{RefQuadrilateral, 1}()
#     qr = QuadratureRule{RefQuadrilateral}(1);
#     qr_facet = FacetQuadratureRule{RefQuadrilateral}(1);
#     kopp_values = ValuesCache(
#         CellValues(qr, ip),
#         FacetValues(qr_facet, ip),
#         InterfaceValues(qr_facet, ip)
#     );
#     @testset "3x3 Quadrilateral corner" begin
#         grid = generate_grid(KoppCell{2, Int}, (3,3))
#         topology = KoppTopology(grid)
#         dh = DofHandler(grid.base_grid)
#         add!(dh, :u, ip)
#         close!(dh);
#         refinement_cache = KoppRefinementCache(grid, topology);
#         sync = LTSAMRSynchronizer(grid, dh, kopp_values, refinement_cache, topology);
#         cells_to_refine = Set([CellIndex(1)]);
#         refine!(grid, topology, refinement_cache, #=kopp_values,=# sync, cells_to_refine)
#         @testset "one level" begin
#             @testset "coords" begin
#                 coords = zeros(Vec{2, Float64}, 4)
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(2))
#                 @test coords ≈ [Vec(-1.0, -1.0), Vec(-4/6, -1.0), Vec(-4/6, -4/6), Vec(-1.0, -4/6)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(3))
#                 @test coords ≈ [Vec(-4/6, -1.0), Vec(-2/6, -1.0), Vec(-2/6, -4/6), Vec(-4/6, -4/6)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(4))
#                 @test coords ≈ [Vec(-4/6, -4/6), Vec(-2/6, -4/6), Vec(-2/6, -2/6), Vec(-4/6, -2/6)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(5))
#                 @test coords ≈ [Vec(-1.0, -4/6), Vec(-4/6, -4/6), Vec(-4/6, -2/6), Vec(-1.0, -2/6)]
#             end
#             @testset "topology" begin
#                 @test getneighborhood(topology, grid, FacetIndex(2,1)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(2,2)) == [FacetIndex(3,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(2,3)) == [FacetIndex(5,1)]
#                 @test getneighborhood(topology, grid, FacetIndex(2,4)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(3,1)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(3,2)) == [FacetIndex(6,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(3,3)) == [FacetIndex(4,1)]
#                 @test getneighborhood(topology, grid, FacetIndex(3,4)) == [FacetIndex(2,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(4,1)) == [FacetIndex(3,3)]
#                 @test getneighborhood(topology, grid, FacetIndex(4,2)) == [FacetIndex(6,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(4,3)) == [FacetIndex(8,1)]
#                 @test getneighborhood(topology, grid, FacetIndex(4,4)) == [FacetIndex(5,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(5,1)) == [FacetIndex(2,3)]
#                 @test getneighborhood(topology, grid, FacetIndex(5,2)) == [FacetIndex(4,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(5,3)) == [FacetIndex(8,1)]
#                 @test getneighborhood(topology, grid, FacetIndex(5,4)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(6,1)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(6,2)) == [FacetIndex(7,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(6,3)) == [FacetIndex(9,1)]
#                 @test getneighborhood(topology, grid, FacetIndex(6,4)) == [FacetIndex(3,2), FacetIndex(4,2)]
#             end
#             @testset "two levels" begin
#                 cells_to_refine = Set([CellIndex(4)]);
#                 refine!(grid, topology, refinement_cache, sync, cells_to_refine)
#                 @testset "coords" begin
#                     coords = zeros(Vec{2, Float64}, 4)
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(5))
#                     @test coords ≈ [Vec(-4/6, -4/6), Vec(-3/6, -4/6), Vec(-3/6, -3/6), Vec(-4/6, -3/6)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(6))
#                     @test coords ≈ [Vec(-3/6, -4/6), Vec(-2/6, -4/6), Vec(-2/6, -3/6), Vec(-3/6, -3/6)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(7))
#                     @test coords ≈ [Vec(-3/6, -3/6), Vec(-2/6, -3/6), Vec(-2/6, -2/6), Vec(-3/6, -2/6)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(8))
#                     @test coords ≈ [Vec(-4/6, -3/6), Vec(-3/6, -3/6), Vec(-3/6, -2/6), Vec(-4/6, -2/6)]
#                 end
#                 @testset "topology" begin
#                     @test getneighborhood(topology, grid, FacetIndex(2,1)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(2,2)) == [FacetIndex(3,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(2,3)) == [FacetIndex(9,1)]
#                     @test getneighborhood(topology, grid, FacetIndex(2,4)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(3,1)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(3,2)) == [FacetIndex(10,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(3,3)) == [FacetIndex(5,1), FacetIndex(6,1)]
#                     @test getneighborhood(topology, grid, FacetIndex(3,4)) == [FacetIndex(2,2)]

#                     @test getneighborhood(topology, grid, FacetIndex(5,1)) == [FacetIndex(3,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(5,2)) == [FacetIndex(6,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(5,3)) == [FacetIndex(8,1)]
#                     @test getneighborhood(topology, grid, FacetIndex(5,4)) == [FacetIndex(9,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(6,1)) == [FacetIndex(3,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(6,2)) == [FacetIndex(10,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(6,3)) == [FacetIndex(7,1)]
#                     @test getneighborhood(topology, grid, FacetIndex(6,4)) == [FacetIndex(5,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(7,1)) == [FacetIndex(6,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(7,2)) == [FacetIndex(10,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(7,3)) == [FacetIndex(12,1)]
#                     @test getneighborhood(topology, grid, FacetIndex(7,4)) == [FacetIndex(8,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(8,1)) == [FacetIndex(5,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(8,2)) == [FacetIndex(7,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(8,3)) == [FacetIndex(12,1)]
#                     @test getneighborhood(topology, grid, FacetIndex(8,4)) == [FacetIndex(9,2)]

#                     @test getneighborhood(topology, grid, FacetIndex(9,1)) == [FacetIndex(2,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(9,2)) == [FacetIndex(8,4), FacetIndex(5,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(9,3)) == [FacetIndex(12,1)]
#                     @test getneighborhood(topology, grid, FacetIndex(9,4)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(10,1)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(10,2)) == [FacetIndex(11,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(10,3)) == [FacetIndex(13,1)]
#                     @test getneighborhood(topology, grid, FacetIndex(10,4)) == [FacetIndex(3,2), FacetIndex(6,2), FacetIndex(7,2)]
#                 end
#                 @testset "coarsening" begin
#                     cells_to_coarsen = Set([CellIndex(4)]);
#                     coarsen!(grid, topology, refinement_cache, sync, cells_to_coarsen, dh);
#                     @testset "coords" begin
#                         coords = zeros(Vec{2, Float64}, 4)
#                         Ferrite.getcoordinates!(coords, grid, CellIndex(4))
#                         @test coords ≈ [Vec(-4/6, -4/6), Vec(-2/6, -4/6), Vec(-2/6, -2/6), Vec(-4/6, -2/6)]
#                         Ferrite.getcoordinates!(coords, grid, CellIndex(5))
#                         @test coords ≈ [Vec(-1.0, -4/6), Vec(-4/6, -4/6), Vec(-4/6, -2/6), Vec(-1.0, -2/6)]
#                     end
#                     @testset "topology" begin
#                         @test getneighborhood(topology, grid, FacetIndex(2,1)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(2,2)) == [FacetIndex(3,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(2,3)) == [FacetIndex(5,1)]
#                         @test getneighborhood(topology, grid, FacetIndex(2,4)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(3,1)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(3,2)) == [FacetIndex(6,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(3,3)) == [FacetIndex(4,1)]
#                         @test getneighborhood(topology, grid, FacetIndex(3,4)) == [FacetIndex(2,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(4,1)) == [FacetIndex(3,3)]
#                         @test getneighborhood(topology, grid, FacetIndex(4,2)) == [FacetIndex(6,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(4,3)) == [FacetIndex(8,1)]
#                         @test getneighborhood(topology, grid, FacetIndex(4,4)) == [FacetIndex(5,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(5,1)) == [FacetIndex(2,3)]
#                         @test getneighborhood(topology, grid, FacetIndex(5,2)) == [FacetIndex(4,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(5,3)) == [FacetIndex(8,1)]
#                         @test getneighborhood(topology, grid, FacetIndex(5,4)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(6,1)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(6,2)) == [FacetIndex(7,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(6,3)) == [FacetIndex(9,1)]
#                         @test getneighborhood(topology, grid, FacetIndex(6,4)) == [FacetIndex(3,2), FacetIndex(4,2)]
#                     end
#                 end
#             end
#         end
#     end
# end

# @testset "Kopp Hexahedron Refinement" begin
#     ip = DiscontinuousLagrange{RefHexahedron, 1}()
#     qr = QuadratureRule{RefHexahedron}(1);
#     qr_facet = FacetQuadratureRule{RefHexahedron}(1);
#     kopp_values = ValuesCache(
#         CellValues(qr, ip),
#         FacetValues(qr_facet, ip),
#         InterfaceValues(qr_facet, ip)
#     );
#     @testset "2x2 Hexahedron corner" begin
#         grid = generate_grid(KoppCell{3, Int}, (2,2,2))
#         topology = KoppTopology(grid)
#         dh = DofHandler(grid.base_grid)
#         add!(dh, :u, ip)
#         close!(dh);
#         refinement_cache = KoppRefinementCache(grid, topology);
#         sync = LTSAMRSynchronizer(grid, dh, kopp_values, refinement_cache, topology);
#         cells_to_refine = Set([CellIndex(1)]);
#         refine!(grid, topology, refinement_cache, sync, cells_to_refine)
#         @testset "one level" begin
#             @testset "coords" begin
#                 coords = zeros(Vec{3, Float64}, 8)
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(2))
#                 @test coords ≈ [Vec(-1.0, -1.0, -1.0), Vec(-0.5, -1.0, -1.0),
#                                 Vec(-0.5, -0.5, -1.0), Vec(-1.0, -0.5, -1.0),
#                                 Vec(-1.0, -1.0, -0.5), Vec(-0.5, -1.0, -0.5),
#                                 Vec(-0.5, -0.5, -0.5), Vec(-1.0, -0.5, -0.5)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(3))
#                 @test coords ≈ [Vec(-0.5, -1.0, -1.0), Vec(0.0, -1.0, -1.0),
#                                 Vec(0.0, -0.5, -1.0), Vec(-0.5, -0.5, -1.0),
#                                 Vec(-0.5, -1.0, -0.5), Vec(0.0, -1.0, -0.5),
#                                 Vec(0.0, -0.5, -0.5), Vec(-0.5, -0.5, -0.5)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(4))
#                 @test coords ≈ [Vec(-0.5, -0.5, -1.0), Vec(0.0, -0.5, -1.0),
#                                 Vec(0.0, 0.0, -1.0), Vec(-0.5, 0.0, -1.0),
#                                 Vec(-0.5, -0.5, -0.5), Vec(0.0, -0.5, -0.5),
#                                 Vec(0.0, 0.0, -0.5), Vec(-0.5, 0.0, -0.5)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(5))
#                 @test coords ≈ [Vec(-1.0, -0.5, -1.0), Vec(-0.5, -0.5, -1.0),
#                                 Vec(-0.5, 0.0, -1.0), Vec(-1.0, 0.0, -1.0),
#                                 Vec(-1.0, -0.5, -0.5), Vec(-0.5, -0.5, -0.5),
#                                 Vec(-0.5, 0.0, -0.5), Vec(-1.0, 0.0, -0.5)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(6))
#                 @test coords ≈ [Vec(-1.0, -1.0, -0.5), Vec(-0.5, -1.0, -0.5),
#                                 Vec(-0.5, -0.5, -0.5), Vec(-1.0, -0.5, -0.5),
#                                 Vec(-1.0, -1.0, -0.0), Vec(-0.5, -1.0, -0.0),
#                                 Vec(-0.5, -0.5, -0.0), Vec(-1.0, -0.5, -0.0)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(7))
#                 @test coords ≈ [Vec(-0.5, -1.0, -0.5), Vec(0.0, -1.0, -0.5),
#                                 Vec(0.0, -0.5, -0.5), Vec(-0.5, -0.5, -0.5),
#                                 Vec(-0.5, -1.0, -0.0), Vec(0.0, -1.0, -0.0),
#                                 Vec(0.0, -0.5, -0.0), Vec(-0.5, -0.5, -0.0)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(8))
#                 @test coords ≈ [Vec(-0.5, -0.5, -0.5), Vec(0.0, -0.5, -0.5),
#                                 Vec(0.0, 0.0, -0.5), Vec(-0.5, 0.0, -0.5),
#                                 Vec(-0.5, -0.5, -0.0), Vec(0.0, -0.5, -0.0),
#                                 Vec(0.0, 0.0, -0.0), Vec(-0.5, 0.0, -0.0)]
#                 Ferrite.getcoordinates!(coords, grid, CellIndex(9))
#                 @test coords ≈ [Vec(-1.0, -0.5, -0.5), Vec(-0.5, -0.5, -0.5),
#                                 Vec(-0.5, 0.0, -0.5), Vec(-1.0, 0.0, -0.5),
#                                 Vec(-1.0, -0.5, -0.0), Vec(-0.5, -0.5, -0.0),
#                                 Vec(-0.5, 0.0, -0.0), Vec(-1.0, 0.0, -0.0)]
#             end
#             @testset "topology" begin
#                 @test getneighborhood(topology, grid, FacetIndex(2,1)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(2,2)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(2,3)) == [FacetIndex(3,5)]
#                 @test getneighborhood(topology, grid, FacetIndex(2,4)) == [FacetIndex(5,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(2,5)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(2,6)) == [FacetIndex(6,1)]

#                 @test getneighborhood(topology, grid, FacetIndex(3,1)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(3,2)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(3,3)) == [FacetIndex(10,5)]
#                 @test getneighborhood(topology, grid, FacetIndex(3,4)) == [FacetIndex(4,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(3,5)) == [FacetIndex(2,3)]
#                 @test getneighborhood(topology, grid, FacetIndex(3,6)) == [FacetIndex(7,1)]

#                 @test getneighborhood(topology, grid, FacetIndex(4,1)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(4,2)) == [FacetIndex(3,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(4,3)) == [FacetIndex(10,5)]
#                 @test getneighborhood(topology, grid, FacetIndex(4,4)) == [FacetIndex(11,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(4,5)) == [FacetIndex(5,3)]
#                 @test getneighborhood(topology, grid, FacetIndex(4,6)) == [FacetIndex(8,1)]

#                 @test getneighborhood(topology, grid, FacetIndex(5,1)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(5,2)) == [FacetIndex(2,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(5,3)) == [FacetIndex(4,5)]
#                 @test getneighborhood(topology, grid, FacetIndex(5,4)) == [FacetIndex(11,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(5,5)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(5,6)) == [FacetIndex(9,1)]

#                 @test getneighborhood(topology, grid, FacetIndex(6,1)) == [FacetIndex(2,6)]
#                 @test getneighborhood(topology, grid, FacetIndex(6,2)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(6,3)) == [FacetIndex(7,5)]
#                 @test getneighborhood(topology, grid, FacetIndex(6,4)) == [FacetIndex(9,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(6,5)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(6,6)) == [FacetIndex(13,1)]

#                 @test getneighborhood(topology, grid, FacetIndex(7,1)) == [FacetIndex(3,6)]
#                 @test getneighborhood(topology, grid, FacetIndex(7,2)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(7,3)) == [FacetIndex(10,5)]
#                 @test getneighborhood(topology, grid, FacetIndex(7,4)) == [FacetIndex(8,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(7,5)) == [FacetIndex(6,3)]
#                 @test getneighborhood(topology, grid, FacetIndex(7,6)) == [FacetIndex(13,1)]

#                 @test getneighborhood(topology, grid, FacetIndex(8,1)) == [FacetIndex(4,6)]
#                 @test getneighborhood(topology, grid, FacetIndex(8,2)) == [FacetIndex(7,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(8,3)) == [FacetIndex(10,5)]
#                 @test getneighborhood(topology, grid, FacetIndex(8,4)) == [FacetIndex(11,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(8,5)) == [FacetIndex(9,3)]
#                 @test getneighborhood(topology, grid, FacetIndex(8,6)) == [FacetIndex(13,1)]

#                 @test getneighborhood(topology, grid, FacetIndex(9,1)) == [FacetIndex(5,6)]
#                 @test getneighborhood(topology, grid, FacetIndex(9,2)) == [FacetIndex(6,4)]
#                 @test getneighborhood(topology, grid, FacetIndex(9,3)) == [FacetIndex(8,5)]
#                 @test getneighborhood(topology, grid, FacetIndex(9,4)) == [FacetIndex(11,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(9,5)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(9,6)) == [FacetIndex(13,1)]

#                 @test getneighborhood(topology, grid, FacetIndex(10,1)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(10,2)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(10,3)) == FacetIndex[]
#                 @test getneighborhood(topology, grid, FacetIndex(10,4)) == [FacetIndex(12,2)]
#                 @test getneighborhood(topology, grid, FacetIndex(10,5)) == [FacetIndex(3,3), FacetIndex(4,3), FacetIndex(8,3), FacetIndex(7,3)]
#                 @test getneighborhood(topology, grid, FacetIndex(10,6)) == [FacetIndex(14,1)]
#             end
#             @testset "two levels" begin
#                 cells_to_refine = Set([CellIndex(8)]);
#                 refine!(grid, topology, refinement_cache, sync, cells_to_refine)
#                 @testset "coords" begin
#                     coords = zeros(Vec{3, Float64}, 8)
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(9))
#                     @test coords ≈ [Vec(-0.5, -0.5, -0.5), Vec(-0.25, -0.5, -0.5),
#                                     Vec(-0.25, -0.25, -0.5), Vec(-0.5, -0.25, -0.5),
#                                     Vec(-0.5, -0.5, -0.25), Vec(-0.25, -0.5, -0.25),
#                                     Vec(-0.25, -0.25, -0.25), Vec(-0.5, -0.25, -0.25)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(10))
#                     @test coords ≈ [Vec(-0.25, -0.5, -0.5), Vec(0.0, -0.5, -0.5),
#                                     Vec(0.0, -0.25, -0.5), Vec(-0.25, -0.25, -0.5),
#                                     Vec(-0.25, -0.5, -0.25), Vec(0.0, -0.5, -0.25),
#                                     Vec(0.0, -0.25, -0.25), Vec(-0.25, -0.25, -0.25)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(11))
#                     @test coords ≈ [Vec(-0.25, -0.25, -0.5), Vec(0.0, -0.25, -0.5),
#                                     Vec(0.0, 0.0, -0.5), Vec(-0.25, 0.0, -0.5),
#                                     Vec(-0.25, -0.25, -0.25), Vec(0.0, -0.25, -0.25),
#                                     Vec(0.0, 0.0, -0.25), Vec(-0.25, 0.0, -0.25)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(12))
#                     @test coords ≈ [Vec(-0.5, -0.25, -0.5), Vec(-0.25, -0.25, -0.5),
#                                     Vec(-0.25, 0.0, -0.5), Vec(-0.5, 0.0, -0.5),
#                                     Vec(-0.5, -0.25, -0.25), Vec(-0.25, -0.25, -0.25),
#                                     Vec(-0.25, 0.0, -0.25), Vec(-0.5, 0.0, -0.25)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(13))
#                     @test coords ≈ [Vec(-0.5, -0.5, -0.25), Vec(-0.25, -0.5, -0.25),
#                                     Vec(-0.25, -0.25, -0.25), Vec(-0.5, -0.25, -0.25),
#                                     Vec(-0.5, -0.5, -0.0), Vec(-0.25, -0.5, -0.0),
#                                     Vec(-0.25, -0.25, -0.0), Vec(-0.5, -0.25, -0.0)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(14))
#                     @test coords ≈ [Vec(-0.25, -0.5, -0.25), Vec(0.0, -0.5, -0.25),
#                                     Vec(0.0, -0.25, -0.25), Vec(-0.25, -0.25, -0.25),
#                                     Vec(-0.25, -0.5, -0.0), Vec(0.0, -0.5, -0.0),
#                                     Vec(0.0, -0.25, -0.0), Vec(-0.25, -0.25, -0.0)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(15))
#                     @test coords ≈ [Vec(-0.25, -0.25, -0.25), Vec(0.0, -0.25, -0.25),
#                                     Vec(0.0, 0.0, -0.25), Vec(-0.25, 0.0, -0.25),
#                                     Vec(-0.25, -0.25, -0.0), Vec(0.0, -0.25, -0.0),
#                                     Vec(0.0, 0.0, -0.0), Vec(-0.25, 0.0, -0.0)]
#                     Ferrite.getcoordinates!(coords, grid, CellIndex(16))
#                     @test coords ≈ [Vec(-0.5, -0.25, -0.25), Vec(-0.25, -0.25, -0.25),
#                                     Vec(-0.25, 0.0, -0.25), Vec(-0.5, 0.0, -0.25),
#                                     Vec(-0.5, -0.25, -0.0), Vec(-0.25, -0.25, -0.0),
#                                     Vec(-0.25, 0.0, -0.0), Vec(-0.5, 0.0, -0.0)]
#                 end
#                 @testset "topology" begin
#                     @test getneighborhood(topology, grid, FacetIndex(2,1)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(2,2)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(2,3)) == [FacetIndex(3,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(2,4)) == [FacetIndex(5,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(2,5)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(2,6)) == [FacetIndex(6,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(3,1)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(3,2)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(3,3)) == [FacetIndex(18,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(3,4)) == [FacetIndex(4,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(3,5)) == [FacetIndex(2,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(3,6)) == [FacetIndex(7,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(4,1)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(4,2)) == [FacetIndex(3,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(4,3)) == [FacetIndex(18,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(4,4)) == [FacetIndex(19,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(4,5)) == [FacetIndex(5,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(4,6)) == [FacetIndex(9,1), FacetIndex(12,1), FacetIndex(11,1), FacetIndex(10,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(5,1)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(5,2)) == [FacetIndex(2,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(5,3)) == [FacetIndex(4,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(5,4)) == [FacetIndex(19,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(5,5)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(5,6)) == [FacetIndex(17,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(6,1)) == [FacetIndex(2,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(6,2)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(6,3)) == [FacetIndex(7,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(6,4)) == [FacetIndex(17,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(6,5)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(6,6)) == [FacetIndex(21,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(7,1)) == [FacetIndex(3,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(7,2)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(7,3)) == [FacetIndex(18,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(7,4)) == [FacetIndex(9,2), FacetIndex(10,2), FacetIndex(14,2), FacetIndex(13,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(7,5)) == [FacetIndex(6,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(7,6)) == [FacetIndex(21,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(9,1)) == [FacetIndex(4,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(9,2)) == [FacetIndex(7,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(9,3)) == [FacetIndex(10,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(9,4)) == [FacetIndex(12,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(9,5)) == [FacetIndex(17,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(9,6)) == [FacetIndex(13,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(10,1)) == [FacetIndex(4,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(10,2)) == [FacetIndex(7,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(10,3)) == [FacetIndex(18,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(10,4)) == [FacetIndex(11,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(10,5)) == [FacetIndex(9,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(10,6)) == [FacetIndex(14,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(11,1)) == [FacetIndex(4,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(11,2)) == [FacetIndex(10,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(11,3)) == [FacetIndex(18,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(11,4)) == [FacetIndex(19,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(11,5)) == [FacetIndex(12,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(11,6)) == [FacetIndex(15,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(12,1)) == [FacetIndex(4,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(12,2)) == [FacetIndex(9,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(12,3)) == [FacetIndex(11,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(12,4)) == [FacetIndex(19,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(12,5)) == [FacetIndex(17,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(12,6)) == [FacetIndex(16,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(13,1)) == [FacetIndex(9,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(13,2)) == [FacetIndex(7,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(13,3)) == [FacetIndex(14,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(13,4)) == [FacetIndex(16,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(13,5)) == [FacetIndex(17,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(13,6)) == [FacetIndex(21,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(14,1)) == [FacetIndex(10,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(14,2)) == [FacetIndex(7,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(14,3)) == [FacetIndex(18,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(14,4)) == [FacetIndex(15,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(14,5)) == [FacetIndex(13,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(14,6)) == [FacetIndex(21,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(15,1)) == [FacetIndex(11,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(15,2)) == [FacetIndex(14,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(15,3)) == [FacetIndex(18,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(15,4)) == [FacetIndex(19,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(15,5)) == [FacetIndex(16,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(15,6)) == [FacetIndex(21,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(16,1)) == [FacetIndex(12,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(16,2)) == [FacetIndex(13,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(16,3)) == [FacetIndex(15,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(16,4)) == [FacetIndex(19,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(16,5)) == [FacetIndex(17,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(16,6)) == [FacetIndex(21,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(17,1)) == [FacetIndex(5,6)]
#                     @test getneighborhood(topology, grid, FacetIndex(17,2)) == [FacetIndex(6,4)]
#                     @test getneighborhood(topology, grid, FacetIndex(17,3)) == [FacetIndex(9,5), FacetIndex(13,5), FacetIndex(16,5), FacetIndex(12,5)]
#                     @test getneighborhood(topology, grid, FacetIndex(17,4)) == [FacetIndex(19,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(17,5)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(17,6)) == [FacetIndex(21,1)]

#                     @test getneighborhood(topology, grid, FacetIndex(18,1)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(18,2)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(18,3)) == FacetIndex[]
#                     @test getneighborhood(topology, grid, FacetIndex(18,4)) == [FacetIndex(20,2)]
#                     @test getneighborhood(topology, grid, FacetIndex(18,5)) == [FacetIndex(3,3), FacetIndex(4,3), FacetIndex(10,3), FacetIndex(11,3), FacetIndex(15,3), FacetIndex(14,3), FacetIndex(7,3)]
#                     @test getneighborhood(topology, grid, FacetIndex(18,6)) == [FacetIndex(22,1)]
#                 end
#                 @testset "coarsening" begin
#                     cells_to_coarsen = Set([CellIndex(8)])
#                     coarsen!(grid, topology, refinement_cache, sync, cells_to_coarsen, dh);
#                     @testset "coords" begin
#                         coords = zeros(Vec{3, Float64}, 8)
#                         Ferrite.getcoordinates!(coords, grid, CellIndex(8))
#                         @test coords ≈ [Vec(-0.5, -0.5, -0.5), Vec(0.0, -0.5, -0.5),
#                                         Vec(0.0, 0.0, -0.5), Vec(-0.5, 0.0, -0.5),
#                                         Vec(-0.5, -0.5, -0.0), Vec(0.0, -0.5, -0.0),
#                                         Vec(0.0, 0.0, -0.0), Vec(-0.5, 0.0, -0.0)]
#                         Ferrite.getcoordinates!(coords, grid, CellIndex(9))
#                         @test coords ≈ [Vec(-1.0, -0.5, -0.5), Vec(-0.5, -0.5, -0.5),
#                                         Vec(-0.5, 0.0, -0.5), Vec(-1.0, 0.0, -0.5),
#                                         Vec(-1.0, -0.5, -0.0), Vec(-0.5, -0.5, -0.0),
#                                         Vec(-0.5, 0.0, -0.0), Vec(-1.0, 0.0, -0.0)]
#                     end
#                     @testset "topology" begin
#                         @test getneighborhood(topology, grid, FacetIndex(2,1)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(2,2)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(2,3)) == [FacetIndex(3,5)]
#                         @test getneighborhood(topology, grid, FacetIndex(2,4)) == [FacetIndex(5,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(2,5)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(2,6)) == [FacetIndex(6,1)]

#                         @test getneighborhood(topology, grid, FacetIndex(3,1)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(3,2)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(3,3)) == [FacetIndex(10,5)]
#                         @test getneighborhood(topology, grid, FacetIndex(3,4)) == [FacetIndex(4,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(3,5)) == [FacetIndex(2,3)]
#                         @test getneighborhood(topology, grid, FacetIndex(3,6)) == [FacetIndex(7,1)]

#                         @test getneighborhood(topology, grid, FacetIndex(4,1)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(4,2)) == [FacetIndex(3,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(4,3)) == [FacetIndex(10,5)]
#                         @test getneighborhood(topology, grid, FacetIndex(4,4)) == [FacetIndex(11,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(4,5)) == [FacetIndex(5,3)]
#                         @test getneighborhood(topology, grid, FacetIndex(4,6)) == [FacetIndex(8,1)]

#                         @test getneighborhood(topology, grid, FacetIndex(5,1)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(5,2)) == [FacetIndex(2,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(5,3)) == [FacetIndex(4,5)]
#                         @test getneighborhood(topology, grid, FacetIndex(5,4)) == [FacetIndex(11,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(5,5)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(5,6)) == [FacetIndex(9,1)]

#                         @test getneighborhood(topology, grid, FacetIndex(6,1)) == [FacetIndex(2,6)]
#                         @test getneighborhood(topology, grid, FacetIndex(6,2)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(6,3)) == [FacetIndex(7,5)]
#                         @test getneighborhood(topology, grid, FacetIndex(6,4)) == [FacetIndex(9,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(6,5)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(6,6)) == [FacetIndex(13,1)]

#                         @test getneighborhood(topology, grid, FacetIndex(7,1)) == [FacetIndex(3,6)]
#                         @test getneighborhood(topology, grid, FacetIndex(7,2)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(7,3)) == [FacetIndex(10,5)]
#                         @test getneighborhood(topology, grid, FacetIndex(7,4)) == [FacetIndex(8,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(7,5)) == [FacetIndex(6,3)]
#                         @test getneighborhood(topology, grid, FacetIndex(7,6)) == [FacetIndex(13,1)]

#                         @test getneighborhood(topology, grid, FacetIndex(8,1)) == [FacetIndex(4,6)]
#                         @test getneighborhood(topology, grid, FacetIndex(8,2)) == [FacetIndex(7,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(8,3)) == [FacetIndex(10,5)]
#                         @test getneighborhood(topology, grid, FacetIndex(8,4)) == [FacetIndex(11,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(8,5)) == [FacetIndex(9,3)]
#                         @test getneighborhood(topology, grid, FacetIndex(8,6)) == [FacetIndex(13,1)]

#                         @test getneighborhood(topology, grid, FacetIndex(9,1)) == [FacetIndex(5,6)]
#                         @test getneighborhood(topology, grid, FacetIndex(9,2)) == [FacetIndex(6,4)]
#                         @test getneighborhood(topology, grid, FacetIndex(9,3)) == [FacetIndex(8,5)]
#                         @test getneighborhood(topology, grid, FacetIndex(9,4)) == [FacetIndex(11,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(9,5)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(9,6)) == [FacetIndex(13,1)]

#                         @test getneighborhood(topology, grid, FacetIndex(10,1)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(10,2)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(10,3)) == FacetIndex[]
#                         @test getneighborhood(topology, grid, FacetIndex(10,4)) == [FacetIndex(12,2)]
#                         @test getneighborhood(topology, grid, FacetIndex(10,5)) == [FacetIndex(8,3), FacetIndex(3,3), FacetIndex(7,3), FacetIndex(4,3)]
#                         @test getneighborhood(topology, grid, FacetIndex(10,6)) == [FacetIndex(14,1)]
#                     end
#                 end
#             end
#         end
#     end
# end

# @testset "Iterators" begin
#     @testset "3x3 Quadrilateral corner" begin
#         grid = generate_grid(KoppCell{2, Int}, (3,3))
#         topology = KoppTopology(grid)
#         ip = DiscontinuousLagrange{RefQuadrilateral, 1}()
#         qr = QuadratureRule{RefQuadrilateral}(1);
#         qr_facet = FacetQuadratureRule{RefQuadrilateral}(1);
#         kopp_values = ValuesCache(
#             CellValues(qr, ip),
#             FacetValues(qr_facet, ip),
#             InterfaceValues(qr_facet, ip)
#         );
#         dh = DofHandler(grid.base_grid)
#         add!(dh, :u, ip)
#         close!(dh);
#         refinement_cache = KoppRefinementCache(grid, topology);
#         sync = LTSAMRSynchronizer(grid, dh, kopp_values, refinement_cache, topology);
#         cells_to_refine = Set([CellIndex(1)]);
#         @testset "0 levels" begin
#             @testset "Topology" begin
#                 topology_vector = [
#                     # ---- cell 1 ---- (Escaped the loop!)
#                     (FacetIndex(1,1), FacetIndex[]),
#                     (FacetIndex(1,2), FacetIndex[FacetIndex(2,4)]),
#                     (FacetIndex(1,3), FacetIndex[FacetIndex(4,1)]),
#                     (FacetIndex(1,4), FacetIndex[]),

#                     # ---- cell 2 ----
#                     (FacetIndex(2,1), FacetIndex[]),
#                     (FacetIndex(2,2), FacetIndex[FacetIndex(3,4)]),
#                     (FacetIndex(2,3), FacetIndex[FacetIndex(5,1)]),
#                     (FacetIndex(2,4), FacetIndex[FacetIndex(1,2)]),

#                     # ---- cell 3 ----
#                     (FacetIndex(3,1), FacetIndex[]),
#                     (FacetIndex(3,2), FacetIndex[]),
#                     (FacetIndex(3,3), FacetIndex[FacetIndex(6,1)]),
#                     (FacetIndex(3,4), FacetIndex[FacetIndex(2,2)]),

#                     # ---- cell 4 ----
#                     (FacetIndex(4,1), FacetIndex[FacetIndex(1,3)]),
#                     (FacetIndex(4,2), FacetIndex[FacetIndex(5,4)]),
#                     (FacetIndex(4,3), FacetIndex[FacetIndex(7,1)]),
#                     (FacetIndex(4,4), FacetIndex[]),

#                     # ---- cell 5 ----
#                     (FacetIndex(5,1), FacetIndex[FacetIndex(2,3)]),
#                     (FacetIndex(5,2), FacetIndex[FacetIndex(6,4)]),
#                     (FacetIndex(5,3), FacetIndex[FacetIndex(8,1)]),
#                     (FacetIndex(5,4), FacetIndex[FacetIndex(4,2)]),

#                     # ---- cell 6 ----
#                     (FacetIndex(6,1), FacetIndex[FacetIndex(3,3)]),
#                     (FacetIndex(6,2), FacetIndex[]),
#                     (FacetIndex(6,3), FacetIndex[FacetIndex(9,1)]),
#                     (FacetIndex(6,4), FacetIndex[FacetIndex(5,2)]),

#                     # ---- cell 7 ----
#                     (FacetIndex(7,1), FacetIndex[FacetIndex(4,3)]),
#                     (FacetIndex(7,2), FacetIndex[FacetIndex(8,4)]),
#                     (FacetIndex(7,3), FacetIndex[]),
#                     (FacetIndex(7,4), FacetIndex[]),

#                     # ---- cell 8 ----
#                     (FacetIndex(8,1), FacetIndex[FacetIndex(5,3)]),
#                     (FacetIndex(8,2), FacetIndex[FacetIndex(9,4)]),
#                     (FacetIndex(8,3), FacetIndex[]),
#                     (FacetIndex(8,4), FacetIndex[FacetIndex(7,2)]),

#                     # ---- cell 9 ----
#                     (FacetIndex(9,1), FacetIndex[FacetIndex(6,3)]),
#                     (FacetIndex(9,2), FacetIndex[]),
#                     (FacetIndex(9,3), FacetIndex[]),
#                     (FacetIndex(9,4), FacetIndex[FacetIndex(8,2)])
#                 ]
#                 i = 1
#                 for (facet, neighborhood) in topology
#                     @test facet == topology_vector[i][1]
#                     @test neighborhood == topology_vector[i][2]
#                     i += 1
#                 end
#             end
#             @testset "CellIterator" begin
#                 cells_vector = [
#                     CellIndex(1)
#                     CellIndex(2)
#                     CellIndex(3)
#                     CellIndex(4)
#                     CellIndex(5)
#                     CellIndex(6)
#                     CellIndex(7)
#                     CellIndex(8)
#                     CellIndex(9)
#                 ]
#                 i = 1
#                 for cc in CellIterator(grid)
#                     @test CellIndex(Ferrite.cellid(cc)) == cells_vector[i]
#                     i+=1
#                 end
#             end
#             @testset "FacetIterator" begin

#             end
#             @testset "InterfaceIterator" begin
#                 interfaces_vector = [
#                     # ---- cell 1 --- (Escaped the loop!)
#                     (FacetIndex(1,2), FacetIndex(2,4)),
#                     (FacetIndex(1,3), FacetIndex(4,1)),

#                     # ---- cell 2 --
#                     (FacetIndex(2,2), FacetIndex(3,4)),
#                     (FacetIndex(2,3), FacetIndex(5,1)),

#                     # ---- cell 3 --
#                     (FacetIndex(3,3), FacetIndex(6,1)),

#                     # ---- cell 4 --
#                     (FacetIndex(4,2), FacetIndex(5,4)),
#                     (FacetIndex(4,3), FacetIndex(7,1)),

#                     # ---- cell 5 --
#                     (FacetIndex(5,2), FacetIndex(6,4)),
#                     (FacetIndex(5,3), FacetIndex(8,1)),

#                     # ---- cell 6 --
#                     (FacetIndex(6,3), FacetIndex(9,1)),

#                     # ---- cell 7 --
#                     (FacetIndex(7,2), FacetIndex(8,4)),
#                     # ---- cell 8 --
#                     (FacetIndex(8,2), FacetIndex(9,4)),
#                 ]
#                 i = 1
#                 for ic in InterfaceIterator(grid, topology)
#                     @test Ferrite.cellid(ic.a) == interfaces_vector[i][1][1]
#                     @test Ferrite.cellid(ic.b) == interfaces_vector[i][2][1]
#                     @test ic.a.current_facet_id == interfaces_vector[i][1][2]
#                     @test ic.b.current_facet_id == interfaces_vector[i][2][2]
#                     i += 1
#                 end
#             end
#         end
#         # @testset "1 levels" begin
#         #     @testset "Topology" begin
#         #         topology_vector = [
#         #             # ---- cell 2 ---- (Escaped the loop!)
#         #             (FacetIndex(1,1), FacetIndex[]),
#         #             (FacetIndex(1,2), FacetIndex[FacetIndex(2,4)]),
#         #             (FacetIndex(1,3), FacetIndex[FacetIndex(4,1)]),
#         #             (FacetIndex(1,4), FacetIndex[]),

#         #             # ---- cell 3 ----
#         #             (FacetIndex(2,1), FacetIndex[]),
#         #             (FacetIndex(2,2), FacetIndex[FacetIndex(3,4)]),
#         #             (FacetIndex(2,3), FacetIndex[FacetIndex(5,1)]),
#         #             (FacetIndex(2,4), FacetIndex[FacetIndex(1,2)]),

#         #             # ---- cell 4 ----
#         #             (FacetIndex(3,1), FacetIndex[]),
#         #             (FacetIndex(3,2), FacetIndex[]),
#         #             (FacetIndex(3,3), FacetIndex[FacetIndex(6,1)]),
#         #             (FacetIndex(3,4), FacetIndex[FacetIndex(2,2)]),

#         #             # ---- cell 5 ----
#         #             (FacetIndex(4,1), FacetIndex[FacetIndex(1,3)]),
#         #             (FacetIndex(4,2), FacetIndex[FacetIndex(5,4)]),
#         #             (FacetIndex(4,3), FacetIndex[FacetIndex(7,1)]),
#         #             (FacetIndex(4,4), FacetIndex[]),

#         #             # ---- cell 6 ----
#         #             (FacetIndex(5,1), FacetIndex[FacetIndex(2,3)]),
#         #             (FacetIndex(5,2), FacetIndex[FacetIndex(6,4)]),
#         #             (FacetIndex(5,3), FacetIndex[FacetIndex(8,1)]),
#         #             (FacetIndex(5,4), FacetIndex[FacetIndex(4,2)]),

#         #             # ---- cell 7 ----
#         #             (FacetIndex(6,1), FacetIndex[FacetIndex(3,3)]),
#         #             (FacetIndex(6,2), FacetIndex[]),
#         #             (FacetIndex(6,3), FacetIndex[FacetIndex(9,1)]),
#         #             (FacetIndex(6,4), FacetIndex[FacetIndex(5,2)]),

#         #             # ---- cell 8 ----
#         #             (FacetIndex(7,1), FacetIndex[FacetIndex(4,3)]),
#         #             (FacetIndex(7,2), FacetIndex[FacetIndex(8,4)]),
#         #             (FacetIndex(7,3), FacetIndex[]),
#         #             (FacetIndex(7,4), FacetIndex[]),

#         #             # ---- cell 9 ----
#         #             (FacetIndex(8,1), FacetIndex[FacetIndex(5,3)]),
#         #             (FacetIndex(8,2), FacetIndex[FacetIndex(9,4)]),
#         #             (FacetIndex(8,3), FacetIndex[]),
#         #             (FacetIndex(8,4), FacetIndex[FacetIndex(7,2)]),

#         #             # ---- cell 10 ----
#         #             (FacetIndex(9,1), FacetIndex[FacetIndex(6,3)]),
#         #             (FacetIndex(9,2), FacetIndex[]),
#         #             (FacetIndex(9,3), FacetIndex[]),
#         #             (FacetIndex(9,4), FacetIndex[FacetIndex(8,2)])
#         #         ]
#         #         i = 1
#         #         for (facet, neighborhood) in topology
#         #             @test facet == topology_vector[i][1]
#         #             @test neighborhood == topology_vector[i][2]
#         #             i += 1
#         #         end
#         #     end
#         #     @testset "CellIterator" begin
#         #         cells_vector = [
#         #             CellIndex(1)
#         #             CellIndex(2)
#         #             CellIndex(3)
#         #             CellIndex(4)
#         #             CellIndex(5)
#         #             CellIndex(6)
#         #             CellIndex(7)
#         #             CellIndex(8)
#         #             CellIndex(9)
#         #         ]
#         #         i = 1
#         #         for cc in CellIterator(grid)
#         #             @test CellIndex(Ferrite.cellid(cc)) == cells_vector[i]
#         #             i+=1
#         #         end
#         #     end
#         #     @testset "FacetIterator" begin

#         #     end
#         #     @testset "InterfaceIterator" begin
#         #         interfaces_vector = [
#         #             # ---- cell 1 --- (Escaped the loop!)
#         #             (FacetIndex(1,2), FacetIndex(2,4)),
#         #             (FacetIndex(1,3), FacetIndex(4,1)),

#         #             # ---- cell 2 --
#         #             (FacetIndex(2,2), FacetIndex(3,4)),
#         #             (FacetIndex(2,3), FacetIndex(5,1)),

#         #             # ---- cell 3 --
#         #             (FacetIndex(3,3), FacetIndex(6,1)),

#         #             # ---- cell 4 --
#         #             (FacetIndex(4,2), FacetIndex(5,4)),
#         #             (FacetIndex(4,3), FacetIndex(7,1)),

#         #             # ---- cell 5 --
#         #             (FacetIndex(5,2), FacetIndex(6,4)),
#         #             (FacetIndex(5,3), FacetIndex(8,1)),

#         #             # ---- cell 6 --
#         #             (FacetIndex(6,3), FacetIndex(9,1)),

#         #             # ---- cell 7 --
#         #             (FacetIndex(7,2), FacetIndex(8,4)),
#         #             # ---- cell 8 --
#         #             (FacetIndex(8,2), FacetIndex(9,4)),
#         #         ]
#         #         i = 1
#         #         for ic in InterfaceIterator(grid, topology)
#         #             @test Ferrite.cellid(ic.a) == interfaces_vector[i][1][1]
#         #             @test Ferrite.cellid(ic.b) == interfaces_vector[i][2][1]
#         #             @test ic.a.current_facet_id == interfaces_vector[i][1][2]
#         #             @test ic.b.current_facet_id == interfaces_vector[i][2][2]
#         #             i += 1
#         #         end
#         #     end
#         # end
#     end
# end

# @testset "2D spiral" begin
#     include("lts_utils.jl")
#     x = 1.025
#     γ = 8.0
#     # refshape = RefTriangle
#     # grid = generate_grid(Triangle, (2,3), Ferrite.Vec((0.0,0.0)), Ferrite.Vec((1.0,1.25)))
#     # transform_coordinates!(grid, x->Ferrite.Vec(
#     #     (1-0.25*(cos(x[2])+0.1)-exp(-3*x[1]),
#     #     1-0.25*(sin(x[1])+0.1)-exp(-4*x[2])
#     # )))

#     refshape = RefQuadrilateral
#     grid = generate_grid(KoppCell{2, Int}, (40,40), Ferrite.Vec((-1.0,0.0)), Ferrite.Vec((1.0,1.25)))
#     # transform_coordinates!(grid.base_grid, x->Ferrite.Vec(
#     #     (1-0.25*(cos(x[2])+0.1)-exp(-3*x[1]),
#     #     1-0.25*(sin(x[1])+0.1)-exp(-4*x[2]),
#     #     # x[3]*x[3]+0.1x[2]
#     #     # x[3]
#     # )))
#     topology = KoppTopology(grid)

#     ip = DiscontinuousLagrange{refshape, 1}()
#     qr = QuadratureRule{refshape}(2)

#     face_qr = FacetQuadratureRule{refshape}(2)
#     cellvalues = CellValues(qr, ip);
#     facetvalues = FacetValues(face_qr, ip)
#     interfacevalues = InterfaceValues(face_qr, ip)

#     dh = DofHandler(grid.base_grid)
#     add!(dh, :u, ip)
#     close!(dh);

#     lts_values = ValuesCache(
#         cellvalues,
#         facetvalues,
#         interfacevalues
#     )

#     refinement_cache = KoppRefinementCache(grid, topology)
#     function spiral_field(x, y; center=(0.,0.5), turns=3.0, clockwise=false)
#         # Calculate offset from center
#         dx = x - center[1]
#         dy = y - center[2]

#         # Convert to polar coordinates
#         r = hypot(dx, dy)  # √(dx² + dy²)
#         θ = atan(dy, dx)   # Angle in [-π, π]

#         # Calculate spiral phase (adjust direction)
#         phase_sign = clockwise ? 1 : -1
#         spiral_phase = θ + phase_sign * 2√2 * turns * π * r

#         return sin(spiral_phase)
#     end
#     function spiral_field(x, y, t; center=(0.,0.5), turns=3.0, clockwise=false, rotation_speed=1.0)
#         # Calculate offset from center
#         dx = x - center[1]
#         dy = y - center[2]

#         # Convert to polar coordinates
#         r = hypot(dx, dy)
#         θ = atan(dy, dx)   # Angle in [-π, π]

#         # Calculate time-dependent rotation (phase shift)
#         time_phase = 2π * rotation_speed * t

#         # Calculate spiral phase with rotation
#         phase_sign = clockwise ? 1 : -1
#         spiral_phase = θ + phase_sign * 2√2 * turns * π * r - time_phase

#         return sin(spiral_phase)
#     end
#     sync = LTSAMRSynchronizer(grid, dh, lts_values, refinement_cache, topology, 0.1)
#     u = sync.data_stores[4].data
#     # Ferrite.apply_analytical!(u, dh, :u, x -> 1/((100000*x[1])^2 + 0.1) + sin(x[2]^3) - 1 )
#     Ferrite.apply_analytical!(u, dh, :u, x -> spiral_field(x[1], x[2]) )
#     Ferrite.apply_analytical!(sync.data_stores_prev[4].data, dh, :u, x -> spiral_field(x[1], x[2]) )
#     # Ferrite.apply_analytical!(sync.data_stores_prev[4].data, dh, :u, x -> 1/((100000*x[1])^2 + 0.1) + sin(x[2]^3) - 1 )
#     needs_refinement = true
#     compute_h(cc::Ferrite.CellCache{<:Any,<:Ferrite.AbstractGrid{1}}) = abs(cc.coords[1][1] - cc.coords[2][1])
#     compute_h(cc::Ferrite.CellCache{<:Any,<:Ferrite.AbstractGrid{2}}) = min(norm(cc.coords[1] - cc.coords[2]), norm(cc.coords[2] - cc.coords[3]), norm(cc.coords[1] - cc.coords[3]))
#     compute_h(cc::Ferrite.CellCache{<:Any,<:Ferrite.AbstractGrid{3}}) = min(
#         norm(cc.coords[1] - cc.coords[2]),
#         norm(cc.coords[1] - cc.coords[3]),
#         norm(cc.coords[1] - cc.coords[4]),
#         norm(cc.coords[2] - cc.coords[3]),
#         norm(cc.coords[2] - cc.coords[4]),
#         norm(cc.coords[3] - cc.coords[4]))
#     compute_h(ic::Ferrite.InterfaceCache) = min(compute_h(ic.a.cc), compute_h(ic.b.cc))

#     refinement_iteration = 0
#     pvd = paraview_collection("amr")

#     for t in 0.0:0.1:1.0
#         needs_refinement = true
#         refinement_iteration = 0
#         while needs_refinement == true
#             needs_refinement = false
#             refinement_set = OrderedSet{CellIndex}()
#             coarsening_set = OrderedSet{CellIndex}()
#             fgrid = to_ferrite_grid(grid)
#             dh2 = deepcopy(sync.dh)
#             resize!(dh2.grid.cells, length(fgrid.cells))
#             dh2.grid.cells .= fgrid.cells
#             resize!(dh2.grid.nodes, length(fgrid.nodes))
#             dh2.grid.nodes .= fgrid.nodes
#             coarsening_vector = zeros(Int, length(fgrid.cells))
#             Ferrite.apply_analytical!(u, dh2, :u, x -> spiral_field(x[1], x[2], t) )
#             # Ferrite.apply_analytical!(sync.data_stores_prev[4].data, dh2, :u, x -> spiral_field(x[1], x[2]) )
#             @time "doing shit" for cc in CellIterator(dh2)
#                 grid.kopp_cells[cellid(cc)].isleaf || continue
#                 h = compute_h(cc)
#                 if mean(u[celldofs(cc)]) * h > 0.005
#                     push!(refinement_set, CellIndex(cellid(cc)))
#                     needs_refinement = true
#                 end
#             end
#             refine!(grid, topology, refinement_cache, sync, refinement_set)


#             fgrid = to_ferrite_grid(grid)
#             dh2 = deepcopy(sync.dh)
#             resize!(dh2.grid.cells, length(fgrid.cells))
#             dh2.grid.cells .= fgrid.cells
#             resize!(dh2.grid.nodes, length(fgrid.nodes))
#             dh2.grid.nodes .= fgrid.nodes
#             coarsening_vector = zeros(Int, length(fgrid.cells))
#             @time "doing shit" for cc in CellIterator(dh2)
#                 grid.kopp_cells[cellid(cc)].isleaf || continue
#                 h = compute_h(cc)
#                 if 0.005 > mean(u[celldofs(cc)]) * h
#                     parent = grid.kopp_cells[cellid(cc)].parent
#                     parent <= 0 && continue
#                     coarsening_vector[parent] += 1
#                 end
#             end
#             @time "doing shit" for cc in CellIterator(dh2)
#                 grid.kopp_cells[cellid(cc)].isleaf && continue
#                 coarsening_vector[cellid(cc)] <4 && continue
#                 push!(coarsening_set, CellIndex(cellid(cc)))
#                 needs_refinement = true
#             end
#             coarsen!(grid, topology, refinement_cache, sync, coarsening_set)
#             refinement_iteration += 1
#             refinement_iteration == 4 && break
#             @warn refinement_iteration
#             # break
#         end
#         fgrid = to_ferrite_grid(grid)
#         dh2 = deepcopy(sync.dh)
#         resize!(dh2.grid.cells, length(fgrid.cells))
#         dh2.grid.cells .= fgrid.cells
#         resize!(dh2.grid.nodes, length(fgrid.nodes))
#         dh2.grid.nodes .= fgrid.nodes
#         VTKGridFile("amr-$t.vtu", fgrid) do vtk
#             write_solution(vtk, dh2, u, "_")
#             pvd[t] = vtk
#         end
#         # VTKGridFile("TTTTTTTT", fgrid) do vtk
#         #     write_solution(vtk, dh2, u, "_")
#         # end;
#     end
#     vtk_save(pvd);

#     # refine!(grid, topology, refinement_cache, sync, Set([
#     #     CellIndex(5),
#     #     CellIndex(14),
#     #     ]))

#     # coarsen!(grid, topology, refinement_cache, sync, Set([
#     #     CellIndex(5),
#     #     # CellIndex(7),
#     #     # CellIndex(8),
#     #     CellIndex(18),
#     #     ]))
#     #     # refine!(grid, topology, refinement_cache, sync, Set([
#     #         # CellIndex(18),
#     #         # ]))
# end

@testset "2D LTS" begin

    include("lts_utils.jl")
    x = 1.025
    γ = 8.0
    # refshape = RefTriangle
    # grid = generate_grid(Triangle, (2,3), Ferrite.Vec((0.0,0.0)), Ferrite.Vec((1.0,1.25)))
    # transform_coordinates!(grid, x->Ferrite.Vec(
    #     (1-0.25*(cos(x[2])+0.1)-exp(-3*x[1]),
    #     1-0.25*(sin(x[1])+0.1)-exp(-4*x[2])
    # )))

    refshape = RefQuadrilateral
    grid = generate_grid(KoppCell{2, Int}, (10,10), Ferrite.Vec((-10.0,-10.0)), Ferrite.Vec((10.0,10.0)))
    # grid = generate_grid(KoppCell{2, Int}, (40,40), Ferrite.Vec((-20.0,0.0)), Ferrite.Vec((20.0,40.25)))
    # grid = generate_grid(KoppCell{2, Int}, (3,3), Ferrite.Vec((-1.0,0.0)), Ferrite.Vec((1.0,1.25)))
    # transform_coordinates!(grid.base_grid, x->Ferrite.Vec(
    #     (1-0.25*(cos(x[2])+0.1)-exp(-3*x[1]),
    #     1-0.25*(sin(x[1])+0.1)-exp(-4*x[2]),
    #     # x[3]*x[3]+0.1x[2]
    #     # x[3]
    # )))
    topology = KoppTopology(grid)

    ip = DiscontinuousLagrange{refshape, 1}()
    qr = QuadratureRule{refshape}(2)

    face_qr = FacetQuadratureRule{refshape}(2)
    cellvalues = CellValues(qr, ip);
    facetvalues = FacetValues(face_qr, ip)
    interfacevalues = InterfaceValues(face_qr, ip)

    dh = DofHandler(grid.base_grid)
    add!(dh, :u, ip)
    close!(dh);

    lts_values = ValuesCache(
        cellvalues,
        facetvalues,
        interfacevalues
    )

    refinement_cache = KoppRefinementCache(grid, topology)
    function spiral_field(x, y; center=(0.,0.), turns=3.0/10, clockwise=false)
        # Calculate offset from center
        dx = x - center[1]
        dy = y - center[2]

        # Convert to polar coordinates
        r = hypot(dx, dy)  # √(dx² + dy²)
        θ = atan(dy, dx)   # Angle in [-π, π]

        # Calculate spiral phase (adjust direction)
        phase_sign = clockwise ? 1 : -1
        spiral_phase = θ + phase_sign * 2√2 * turns * π * r

        return sin(spiral_phase)
    end

    sync = LTSAMRSynchronizer(grid, dh, lts_values, refinement_cache, topology, 0.1)
    u = sync.data_stores[4].data
    # Ferrite.apply_analytical!(u, dh, :u, x -> 1/((100000*x[1])^2 + 0.1) + sin(x[2]^3) - 1 )
    # Ferrite.apply_analytical!(u, dh, :u, x -> 1 - (x[1]^2 + x[2]^2)/200)
    # Ferrite.apply_analytical!(u, dh, :u, x -> spiral_field(x[1], x[2]) )
    Ferrite.apply_analytical!(u, dh, :u, x -> x[1] > 0 && x[2] > 0 ? 10.0 : 0.0 )
    # Ferrite.apply_analytical!(sync.data_stores_prev[4].data, dh, :u, x -> spiral_field(x[1], x[2]) )
    # Ferrite.apply_analytical!(sync.data_stores_prev[4].data, dh, :u, x -> 1/((100000*x[1])^2 + 0.1) + sin(x[2]^3) - 1 )
    needs_refinement = true
    compute_h(cc::Ferrite.CellCache{<:Any,<:Ferrite.AbstractGrid{1}}) = abs(cc.coords[1][1] - cc.coords[2][1])
    compute_h(cc::Ferrite.CellCache{<:Any,<:Ferrite.AbstractGrid{2}}) = min(norm(cc.coords[1] - cc.coords[2]), norm(cc.coords[2] - cc.coords[3]), norm(cc.coords[1] - cc.coords[3]))
    compute_h(cc::Ferrite.CellCache{<:Any,<:Ferrite.AbstractGrid{3}}) = min(
        norm(cc.coords[1] - cc.coords[2]),
        norm(cc.coords[1] - cc.coords[3]),
        norm(cc.coords[1] - cc.coords[4]),
        norm(cc.coords[2] - cc.coords[3]),
        norm(cc.coords[2] - cc.coords[4]),
        norm(cc.coords[3] - cc.coords[4]))
    compute_h(ic::Ferrite.InterfaceCache) = min(compute_h(ic.a.cc), compute_h(ic.b.cc))

    refinement_iteration = 0
    pvd = paraview_collection("amr")

    # ∂ₜu = K u + f
    # -> M(uₙ₊₁ - uₙ)/Δt = K uₙ + f
    # -> uₙ₊₁ - uₙ = Δt*inv(M)(K uₙ + f)
    # -> uₙ₊₁ = uₙ + Δt*inv(M)(K uₙ + f)
    # MinvK = Minv*-K

    # λ,_ = eigen(MinvK)
    # Δt = 1.0/maximum(-λ)
    # @show Δt

    # ndofs_per_element = getnbasefunctions(ip)
    # Δλlocal = [1.0/maximum(-eigvals(MinvK[(ndofs_per_element*n+1):ndofs_per_element*(n+1), (ndofs_per_element*n+1):ndofs_per_element*(n+1)])) for n in 0:(length(grid.cells)-1)]
    # Δσlocal = [1.0/maximum(svdvals(MinvK[(ndofs_per_element*n+1):ndofs_per_element*(n+1), :])) for n in 0:(length(grid.cells)-1)]

    # solve_global(Δt, 1.0)
    # solve_global(minimum(Δλlocal), 1.0)
    # solve_global(minimum(Δσlocal), 1.0)

    function solve_alts(dh, topology, Δte, T)
        grid = dh.grid
        # uₙ   = 100*rand(ndofs(dh))
        uₙ₊₁ = copy(uₙ)
        ncells = length(grid.cells)

        # Assume nodal interpolation
        # xvis = vcat([[Makie.Point2f(grid.nodes[node].x.data) for node in cell.nodes] for cell in grid.cells]...)
        te = zeros(ncells)
        teprev = zeros(ncells)
        while any(te .< T)
            # Which element next
            nei = argmin(te)
            # Local time of the element
            t = te[nei]
            # Set time step length
            Δt = Δte[nei]
            if t + Δt > T
                Δt = T-t
            end

            # Store current solution
            dofs = celldofs(dh, nei)
            uₙ[dofs] .= uₙ₊₁[dofs]
            # TODO be smarter.
            uhelp = zeros(ndofs(dh))
            uhelp[dofs] .= uₙ[dofs]
            # Perform single step
            for lfi ∈ 1:Ferrite.nfaces(grid.cells[nei])
                for face_neighbor ∈ Ferrite.getneighborhood(topology, grid, FaceIndex(nei, lfi))
                    dofsnbr = celldofs(dh, face_neighbor[1])
                    uhelp[dofsnbr] .= interpolate_linear(teprev[face_neighbor[1]], te[face_neighbor[1]], uₙ[dofsnbr], uₙ₊₁[dofsnbr], t)
                end
            end
            uₙ₊₁[dofs] .= uₙ[dofs] + Δt * MinvK[dofs,:]*uhelp

            # Update local time
            teprev[nei] = te[nei]
            te[nei] += Δt

            # Diverged?
            if norm(uₙ₊₁) > 1e4 || any(isnan.(uₙ₊₁))
                @show "Broken at $t with $(norm(uₙ₊₁))"
                break
            end
            # notify(uvis)
            # sleep(0.0001)
            # Monitor if solution grows
            @show nei, te[nei], extrema(uₙ₊₁)
        end
        notify(uvis)
        @show "Done."
    end

    @time "ALL FOR ONE" begin
        # u   .= 100*rand(ndofs(dh))
        # uₙ₊₁ = copy(uₙ)
        # Δt = 0.00013019726292697877
        Δt = 0.00025
        Δt = 0.0025
        # Δt = 0.0000001
        T = 100*Δt
        T = 1.0
        ndofs_cell = 4

        for t ∈ 0.0:Δt:T
            u_new = copy(u)
            error_vector = zeros(length(grid.kopp_cells))
            # fgrid = to_ferrite_grid(grid)
            # dh2 = deepcopy(sync.dh)
            # resize!(dh2.grid.cells, length(fgrid.cells))
            # dh2.grid.cells .= fgrid.cells
            # resize!(dh2.grid.nodes, length(fgrid.nodes))
            # dh2.grid.nodes .= fgrid.nodes
            K = [zeros(Float64, ndofs_cell, (1 + sum(topology.cell_facet_neighbors_length[:, cell])) * ndofs_cell) for cell in 1:length(grid.kopp_cells)]
            dofs_map = [zeros(Int,(1 + sum(topology.cell_facet_neighbors_length[:, cell])) * ndofs_cell) for cell in 1:length(grid.kopp_cells)]
            for cc in CellIterator(grid)
                grid.kopp_cells[cellid(cc)].isleaf || continue
                ndofs = 4
                K[cellid(cc)][1:ndofs,1:ndofs] .= sync.data_stores[2].data[:,:,cellid(cc)]
            end

            interface_index = 1
            for ic in InterfaceIterator(grid, topology)
                cell_a = cellid(ic.a)
                cell_b = cellid(ic.b)
                # dofs_a = celldofs(dh, cell_a)
                # dofs_b = celldofs(dh, cell_b)
                ndofs = 4
                # K[cell_a][1:ndofs,1:ndofs] .= sync.data_stores[2].data[1:ndofs,1:ndofs,cell_a]
                K[cell_a][1:ndofs,1:ndofs] .+= sync.data_stores[3].data[1:ndofs,1:ndofs,interface_index]
                dofs_map[cell_a][1:ndofs] .= celldofs(dh, cell_a)

                neighbor_idx = 1
                found = false
                for facet in 1:4
                    # if found
                    #     break
                    # end
                    neighborhood = getneighborhood(topology, FacetIndex(cell_a, facet))
                    if cell_a == 2
                        # @info facet
                        # @show getneighborhood(topology, FacetIndex(2, facet)) grid.kopp_cells[2] get_refinement_level(grid.kopp_cells[2])
                        # @show getneighborhood(topology, FacetIndex(2185, facet)) grid.kopp_cells[2185] get_refinement_level(grid.kopp_cells[2185])
                        # @show getneighborhood(topology, FacetIndex(2190, facet)) grid.kopp_cells[2190] get_refinement_level(grid.kopp_cells[2190])
                        # @show grid.kopp_cells[2549] get_refinement_level(grid.kopp_cells[2549])
                        # @show grid.kopp_cells[2550] get_refinement_level(grid.kopp_cells[2550])
                    end
                    for (i, neighbor) in enumerate(neighborhood)
                        if neighbor[1] == cell_b
                            dof_start = ndofs + 1 + (neighbor_idx-1)*ndofs
                            dof_end = ndofs + neighbor_idx*ndofs
                            # K[cell_a][1:ndofs,dof_start:dof_end] += sync.data_stores[3].data[ndofs + 1:end, 1:ndofs, interface_index]
                            K[cell_a][1:ndofs,dof_start:dof_end] += sync.data_stores[3].data[1:ndofs, ndofs + 1:end, interface_index]
                            dofs_map[cell_a][dof_start:dof_end] .= celldofs(dh, neighbor[1])
                            if any(celldofs(dh, neighbor[1]) .== 0)
                                @error "NIIIII"
                            end
                            if !(grid.kopp_cells[neighbor[1]].isleaf)
                                @error "NIIIII"
                            end
                            found = true
                            break
                        end
                        neighbor_idx += 1
                    end
                end

                # K[cell_b][1:ndofs,1:ndofs] .= sync.data_stores[2].data[1:end, 1:end,cell_b]
                K[cell_b][1:ndofs,1:ndofs] .+= sync.data_stores[3].data[ndofs+1:end, ndofs+1:end,interface_index]
                # display(sync.data_stores[3].data[:, :,interface_index])
                dofs_map[cell_b][1:ndofs] .= celldofs(dh, cell_b)

                neighbor_idx = 1
                found = false
                for facet in 1:4
                    if found
                        break
                    end
                    neighborhood = getneighborhood(topology, FacetIndex(cell_b, facet))
                    # if cell_b == 2551
                    #     @info facet
                    #     @show getneighborhood(topology, FacetIndex(2542, facet)) grid.kopp_cells[2542] get_refinement_level(grid.kopp_cells[2542])
                    #     @show getneighborhood(topology, FacetIndex(2551, facet)) grid.kopp_cells[2551] get_refinement_level(grid.kopp_cells[2551])
                    #     @show getneighborhood(topology, FacetIndex(2554, facet)) grid.kopp_cells[2554] get_refinement_level(grid.kopp_cells[2554])
                    #     @show grid.kopp_cells[2549] get_refinement_level(grid.kopp_cells[2549])
                    #     @show grid.kopp_cells[2550] get_refinement_level(grid.kopp_cells[2550])
                    # end
                    for (i, neighbor) in enumerate(neighborhood)
                        if neighbor[1] == cell_a
                            dof_start = ndofs + 1 + (neighbor_idx-1)*ndofs
                            dof_end = ndofs + neighbor_idx*ndofs
                            # K[cell_b][1:ndofs,dof_start:dof_end] += sync.data_stores[3].data[1:ndofs, ndofs + 1:end, interface_index]
                            K[cell_b][1:ndofs,dof_start:dof_end] += sync.data_stores[3].data[ndofs + 1:end, 1:ndofs, interface_index]
                            dofs_map[cell_b][dof_start:dof_end] .= celldofs(dh, neighbor[1])
                            if any(celldofs(dh, neighbor[1]) .== 0)
                                @error "NIIIII"
                            end
                            if !(grid.kopp_cells[neighbor[1]].isleaf)
                                @error "NIIIII"
                            end
                            found = true
                            break
                        end
                        neighbor_idx += 1
                    end
                end
                # display(sync.data_stores[3].data[:, :, interface_index])
                interface_index += 1
                reinit!(interfacevalues, ic, topology)
                estimate_kelly_interface!(Float64, error_vector, u[[celldofs(dh, cell_a)..., celldofs(dh, cell_b)...]], ic, sync.values_cache.interface_values)
            end
            for cc in CellIterator(grid)
                grid.kopp_cells[cellid(cc)].isleaf || continue
                # @info cellid(cc) dofs_map[cellid(cc)]
                # display((@view sync.data_stores[1].data[:,:,cellid(cc)]))
                Minv = inv((@view sync.data_stores[1].data[:,:,cellid(cc)]))
                u_new[celldofs(dh, cellid(cc))] .= u[celldofs(dh, cellid(cc))] + Δt * Minv * (-K[cellid(cc)]*u[dofs_map[cellid(cc)]])
                # λ,_ = eigen(Minv * (-K[:,:, cellid(cc)]))
            end
            display(extrema(u_new))
            u .= u_new
            if norm(u) > 1e3
                @show "Broken at $t"
                break
            end
            needs_refinement = true
            refinement_iteration = 0
            threshold = 0.8
            while needs_refinement == true
                needs_refinement = false
                refinement_set = OrderedSet{CellIndex}()
                coarsening_set = OrderedSet{CellIndex}()
                # fgrid = to_ferrite_grid(grid)
                # dh2 = deepcopy(sync.dh)
                # resize!(dh2.grid.cells, length(fgrid.cells))
                # dh2.grid.cells .= fgrid.cells
                # resize!(dh2.grid.nodes, length(fgrid.nodes))
                # dh2.grid.nodes .= fgrid.nodes
                coarsening_vector = zeros(Int, length(grid.kopp_cells))
                @info extrema(error_vector)
                for cc in CellIterator(grid)
                    grid.kopp_cells[cellid(cc)].isleaf || continue
                    h = compute_h(cc)
                    if sqrt(error_vector[cellid(cc)]) > threshold
                        # mean(u[celldofs(cc)]) * h > threshold
                        get_refinement_level(grid.kopp_cells[cellid(cc)]) >= 3 && continue
                        push!(refinement_set, CellIndex(cellid(cc)))
                        needs_refinement = true
                    end
                end
                refine!(grid, topology, refinement_cache, sync, refinement_set)
                resize!(error_vector, length(grid.kopp_cells))
                error_vector .= 0.0
                for ic in InterfaceIterator(grid, topology)
                    cell_a = cellid(ic.a)
                    cell_b = cellid(ic.b)
                    reinit!(interfacevalues, ic, topology)
                    estimate_kelly_interface!(Float64, error_vector, u[[celldofs(dh, cell_a)..., celldofs(dh, cell_b)...]], ic, sync.values_cache.interface_values)
                end

                # fgrid = to_ferrite_grid(grid)
                # dh2 = deepcopy(sync.dh)
                # resize!(dh2.grid.cells, length(fgrid.cells))
                # dh2.grid.cells .= fgrid.cells
                # resize!(dh2.grid.nodes, length(fgrid.nodes))
                # dh2.grid.nodes .= fgrid.nodes
                coarsening_vector = zeros(Int, length(grid.kopp_cells))
                for cc in CellIterator(grid)
                    grid.kopp_cells[cellid(cc)].isleaf || continue
                    # h = compute_h(cc)
                    if threshold/10 > sqrt(error_vector[cellid(cc)])
                        parent = grid.kopp_cells[cellid(cc)].parent
                        parent <= 0 && continue
                        coarsening_vector[parent] += 1
                    end
                end
                for cc in CellIterator(grid)
                    grid.kopp_cells[cellid(cc)].isleaf && continue
                    coarsening_vector[cellid(cc)] <4 && continue
                    push!(coarsening_set, CellIndex(cellid(cc)))
                    needs_refinement = true
                end
                coarsen!(grid, topology, refinement_cache, sync, coarsening_set)
                refinement_iteration += 1
                refinement_iteration == 1 && break
                @warn refinement_iteration
            end
            begin
                fgrid = to_ferrite_grid(grid)
                dh2 = deepcopy(sync.dh)
                resize!(dh2.grid.cells, length(fgrid.cells))
                dh2.grid.cells .= fgrid.cells
                resize!(dh2.grid.nodes, length(fgrid.nodes))
                dh2.grid.nodes .= fgrid.nodes
            end
            resize!(error_vector, length(grid.kopp_cells))
            error_vector .= 0.0
            for ic in InterfaceIterator(grid, topology)
                cell_a = cellid(ic.a)
                cell_b = cellid(ic.b)
                reinit!(interfacevalues, ic, topology)
                estimate_kelly_interface!(Float64, error_vector, u[[celldofs(dh, cell_a)..., celldofs(dh, cell_b)...]], ic, sync.values_cache.interface_values)
            end
            VTKGridFile("amr-$t.vtu", fgrid; write_discontinuous = true) do vtk
                write_solution(vtk, dh2, u, "_")
                write_cell_data(vtk, error_vector, "__")
                pvd[t] = vtk
            end
            # vtk_save(pvd);
                @warn t
            # u .= uₙ₊₁
        end
    end
    threshold = 0.5
    for t in 0.0:0.1:1.0

        # VTKGridFile("TTTTTTTT", fgrid) do vtk
        #     write_solution(vtk, dh2, u, "_")
        # end;
    end
    vtk_save(pvd);

    # refine!(grid, topology, refinement_cache, sync, Set([
    #     CellIndex(5),
    #     CellIndex(14),
    #     ]))

    # coarsen!(grid, topology, refinement_cache, sync, Set([
    #     CellIndex(5),
    #     # CellIndex(7),
    #     # CellIndex(8),
    #     CellIndex(18),
    #     ]))
    #     # refine!(grid, topology, refinement_cache, sync, Set([
    #         # CellIndex(18),
    #         # ]))
end
