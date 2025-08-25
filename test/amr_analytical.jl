using Ferrite
using Tensors
using Test
using Logging
using ForwardDiff
import SHA
using Random
using LinearAlgebra:norm2, svdvals
using SparseArrays
using StaticArrays
using OrderedCollections
using WriteVTK
import Metis
using HCubature: hcubature, hquadrature
using Unroll

include("/home/amohamed/.julia/dev/Ferrite/src/Grid/Adaptivity/kopp.jl")
include("lts_utils.jl")
begin
    x = 1.025
    γ = 8.0
    # refshape = RefTriangle
    # grid = generate_grid(Triangle, (2,3), Ferrite.Vec((0.0,0.0)), Ferrite.Vec((1.0,1.25)))
    # transform_coordinates!(grid, x->Ferrite.Vec(
    #     (1-0.25*(cos(x[2])+0.1)-exp(-3*x[1]),
    #     1-0.25*(sin(x[1])+0.1)-exp(-4*x[2])
    # )))

    refshape = RefQuadrilateral
    grid = generate_grid(KoppCell{2, Int}, (40,40), Ferrite.Vec((-1.0,0.0)), Ferrite.Vec((1.0,1.25)))
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
    function spiral_field(x, y; center=(0.,0.5), turns=3.0, clockwise=false)
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
    function spiral_field(x, y, t; center=(0.,0.5), turns=3.0, clockwise=false, rotation_speed=1.0)
        # Calculate offset from center
        dx = x - center[1]
        dy = y - center[2]

        # Convert to polar coordinates
        r = hypot(dx, dy)
        θ = atan(dy, dx)   # Angle in [-π, π]

        # Calculate time-dependent rotation (phase shift)
        time_phase = 2π * rotation_speed * t

        # Calculate spiral phase with rotation
        phase_sign = clockwise ? 1 : -1
        spiral_phase = θ + phase_sign * 2√2 * turns * π * r - time_phase

        return sin(spiral_phase)
    end
    sync = LTSAMRSynchronizer(grid, dh, lts_values, refinement_cache, topology, 0.1)
    u = sync.data_stores.solution_vector.data
    # Ferrite.apply_analytical!(u, dh, :u, x -> 1/((100000*x[1])^2 + 0.1) + sin(x[2]^3) - 1 )
    Ferrite.apply_analytical!(u, dh, :u, x -> spiral_field(x[1], x[2]) )
    Ferrite.apply_analytical!(sync.data_stores_prev.solution_vector.data, dh, :u, x -> spiral_field(x[1], x[2]) )
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
    pvd = paraview_collection("analytical_spiral")

    for t in 0.0:0.1:1.0
        needs_refinement = true
        refinement_iteration = 0
        while needs_refinement == true
            needs_refinement = false
            refinement_set = OrderedSet{CellIndex}()
            coarsening_set = OrderedSet{CellIndex}()
            sizehint!(refinement_set, length(grid.kopp_cells))
            sizehint!(coarsening_set, length(grid.kopp_cells))
            fgrid = to_ferrite_grid(grid)
            dh2 = deepcopy(sync.dh)
            resize!(dh2.grid.cells, length(fgrid.cells))
            dh2.grid.cells .= fgrid.cells
            resize!(dh2.grid.nodes, length(fgrid.nodes))
            dh2.grid.nodes .= fgrid.nodes
            Ferrite.apply_analytical!(u, dh2, :u, x -> spiral_field(x[1], x[2], t) )
            # Ferrite.apply_analytical!(sync.data_stores_prev[4].data, dh2, :u, x -> spiral_field(x[1], x[2]) )
            @inbounds for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                grid.kopp_cells[cellid(cc)].isleaf || continue
                h = compute_h(cc)
                if mean(u[celldofs(dh, cellid(cc))]) * h > 0.005
                    # get_refinement_level(grid.kopp_cells[cellid(cc)]) >= 2 && continue
                    push!(refinement_set, CellIndex(cellid(cc)))
                    needs_refinement = true
                end
            end
            refine!(grid, topology, refinement_cache, sync, refinement_set)

            coarsening_vector = zeros(Int, length(grid.kopp_cells))
            for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                grid.kopp_cells[cellid(cc)].isleaf || continue
                h = compute_h(cc)
                if -0.005 > mean(u[celldofs(dh, cellid(cc))]) / h
                    parent = grid.kopp_cells[cellid(cc)].parent
                    parent <= 0 && continue
                    coarsening_vector[parent] += 1
                end
            end
            for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                grid.kopp_cells[cellid(cc)].isleaf && continue
                coarsening_vector[cellid(cc)] <4 && continue
                push!(coarsening_set, CellIndex(cellid(cc)))
                needs_refinement = true
            end
            coarsen!(grid, topology, refinement_cache, sync, coarsening_set)
            refinement_iteration += 1
            refinement_iteration == 4 && break
            @warn refinement_iteration
            # break
        end
        fgrid = to_ferrite_grid(grid)
        dh2 = deepcopy(sync.dh)
        resize!(dh2.grid.cells, length(fgrid.cells))
        dh2.grid.cells .= fgrid.cells
        resize!(dh2.grid.nodes, length(fgrid.nodes))
        dh2.grid.nodes .= fgrid.nodes
        VTKGridFile("analytical_spiral/amr-$t.vtu", fgrid) do vtk
            write_solution(vtk, dh2, u, "_")
            pvd[t] = vtk
        end
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
