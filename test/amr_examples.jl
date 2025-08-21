using Ferrite
using Tensors
using Test
using Logging
using ForwardDiff
import SHA
using Random
using LinearAlgebra:norm2
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
    grid = generate_grid(KoppCell{2, Int}, (40,40), Ferrite.Vec((-10.0,-10.0)), Ferrite.Vec((10.0,10.0)))
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
    function spiral_field(x, y; center=(0.,0.), turns=3.0/40, clockwise=false)
        # Calculate offset from center
        dx = x - center[1]
        dy = y - center[2]

        # Convert to polar coordinates
        r = hypot(dx, dy)  # √(dx² + dy²)
        θ = atan(dy, dx)   # Angle in [-π, π]

        # Calculate spiral phase (adjust direction)
        phase_sign = clockwise ? 1 : -1
        spiral_phase = θ + phase_sign * 2√2 * turns * π * r

        return sin(spiral_phase) >0 ? 1. : 0.
    end

    sync = LTSAMRSynchronizer(grid, dh, lts_values, refinement_cache, topology, 0.1)
    u = sync.data_stores[4].data
    # Ferrite.apply_analytical!(u, dh, :u, x -> 1/((100000*x[1])^2 + 0.1) + sin(x[2]^3) - 1 )
    # Ferrite.apply_analytical!(u, dh, :u, x -> 1 - (x[1]^2 + x[2]^2)/200)
    Ferrite.apply_analytical!(u, dh, :u, x -> spiral_field(x[1], x[2]) )
    # Ferrite.apply_analytical!(u, dh, :u, x -> x[1] > 0 && x[2] > 0 ? 10.0 : 0.0 )
    # Ferrite.apply_analytical!(sync.data_stores_prev[4].data, dh, :u, x -> spiral_field(x[1], x[2]) )
    # Ferrite.apply_analytical!(sync.data_stores_prev[4].data, dh, :u, x -> 1/((100000*x[1])^2 + 0.1) + sin(x[2]^3) - 1 )
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
    function solve_direct(grid, topology, dh, refinement_cache, sync, dt, io_enabled::Bool = false)
        pvd = io_enabled ? paraview_collection("amr") : nothing


        # u   .= 100*rand(ndofs(dh))
        # uₙ₊₁ = copy(uₙ)
        # Δt = 0.00013019726292697877
        Δt = 0.00025
        Δt = 0.0025
        # Δt = 0.0000001
        T = 100*Δt
        T = 1.0
        ndofs_cell = 4
        u = sync.data_stores[4].data
        u_new = copy(u)
        error_vector = zeros(length(grid.kopp_cells))

        @views for t ∈ 0.0:Δt:T
            u_new .= 0.0
            error_vector .= 0.0
            K = [zeros(Float64, ndofs_cell, (1 + sum(topology.cell_facet_neighbors_length[:, cell])) * ndofs_cell) for cell in 1:length(grid.kopp_cells)]
            dofs_map = [zeros(Int,(1 + sum(topology.cell_facet_neighbors_length[:, cell])) * ndofs_cell) for cell in 1:length(grid.kopp_cells)]
            @inbounds @views for cc in CellIterator(grid)
                grid.kopp_cells[cellid(cc)].isleaf || continue
                ndofs = 4
                K[cellid(cc)][1:ndofs,1:ndofs] .= sync.data_stores[2].data[:,:,cellid(cc)]
            end


            @views for cc in CellIterator(grid)
                grid.kopp_cells[cellid(cc)].isleaf || continue
                Minv = inv((sync.data_stores[1].data[:,:,cellid(cc)]))

                celldofs!(dofs_temp_storage, dh, cellid(cc))
                u_new[dofs_temp_storage] .= u[dofs_temp_storage] + Δt * Minv * (-K[cellid(cc)]*u[dofs_map[cellid(cc)]])
            end
            display(extrema(u_new))
            u .= u_new
            if norm(u) > 1e3
                @show "Broken at $t"
                break
            end
            needs_refinement = true
            refinement_iteration = 0
            threshold = 0.2
            while needs_refinement == true
                needs_refinement = false
                refinement_set = OrderedSet{CellIndex}()
                coarsening_set = OrderedSet{CellIndex}()
                sizehint!(refinement_set, length(grid.kopp_cells))
                sizehint!(coarsening_set, length(grid.kopp_cells))
                coarsening_vector = zeros(Int, length(grid.kopp_cells))
                @info extrema(error_vector)
                @inbounds for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                    grid.kopp_cells[cellid(cc)].isleaf || continue
                    if sqrt(error_vector[cellid(cc)]) > threshold
                        get_refinement_level(grid.kopp_cells[cellid(cc)]) >= 1 && continue
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
                    # TODO: urgent why reinit allocates
                    reinit!(interfacevalues, ic, topology)
                    celldofs!(( dofs_temp_storage2[1:ndofs_cell]), dh, cell_a)
                    celldofs!(( dofs_temp_storage2[ndofs_cell + 1:end]), dh, cell_b)

                    estimate_kelly_interface!(Float64, error_vector, u[dofs_temp_storage2], ic, sync.values_cache.interface_values)
                end
                coarsening_vector = zeros(Int, length(grid.kopp_cells))
                for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                    grid.kopp_cells[cellid(cc)].isleaf || continue
                    if threshold/10 > sqrt(error_vector[cellid(cc)])
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
                refinement_iteration == 1 && break
                @warn refinement_iteration
            end
            resize!(error_vector, length(grid.kopp_cells))
            resize!(u_new, length(u))
            error_vector .= 0.0
            for ic in InterfaceIterator(grid, topology)
                cell_a = cellid(ic.a)
                cell_b = cellid(ic.b)
                reinit!(interfacevalues, ic, topology)
                celldofs!(( dofs_temp_storage2[1:ndofs_cell]), dh, cell_a)
                celldofs!(( dofs_temp_storage2[ndofs_cell + 1:end]), dh, cell_b)

                estimate_kelly_interface!(Float64, error_vector, u[dofs_temp_storage2], ic, sync.values_cache.interface_values)
            end
            if io_enabled
                fgrid = to_ferrite_grid(grid)
                dh2 = deepcopy(sync.dh)
                resize!(dh2.grid.cells, length(fgrid.cells))
                dh2.grid.cells .= fgrid.cells
                resize!(dh2.grid.nodes, length(fgrid.nodes))
                dh2.grid.nodes .= fgrid.nodes

                VTKGridFile("amr-$t.vtu", fgrid; write_discontinuous = true) do vtk
                    write_solution(vtk, dh2, u, "_")
                    write_cell_data(vtk, error_vector, "__")
                    pvd[t] = vtk
                end

            end
        end


        # Assume nodal interpolation
        # xvis = vcat([[Makie.Point2f(grid.nodes[node].x.data) for node in cell.nodes] for cell in grid.cells]...)
        # te = zeros(ncells)
        # teprev = zeros(ncells)
        # while any(te .< T)
        #     # Which element next
        #     nei = argmin(te)
        #     # Local time of the element
        #     t = te[nei]
        #     # Set time step length
        #     Δt = Δte[nei]
        #     if t + Δt > T
        #         Δt = T-t
        #     end

        #     # Store current solution
        #     dofs = celldofs(dh, nei)
        #     uₙ[dofs] .= uₙ₊₁[dofs]
        #     # TODO be smarter.
        #     uhelp = zeros(ndofs(dh))
        #     uhelp[dofs] .= uₙ[dofs]
        #     # Perform single step
        #     for lfi ∈ 1:Ferrite.nfaces(grid.cells[nei])
        #         for face_neighbor ∈ Ferrite.getneighborhood(topology, grid, FaceIndex(nei, lfi))
        #             dofsnbr = celldofs(dh, face_neighbor[1])
        #             uhelp[dofsnbr] .= interpolate_linear(teprev[face_neighbor[1]], te[face_neighbor[1]], uₙ[dofsnbr], uₙ₊₁[dofsnbr], t)
        #         end
        #     end
        #     uₙ₊₁[dofs] .= uₙ[dofs] + Δt * MinvK[dofs,:]*uhelp

        #     # Update local time
        #     teprev[nei] = te[nei]
        #     te[nei] += Δt

        #     # Diverged?
        #     if norm(uₙ₊₁) > 1e4 || any(isnan.(uₙ₊₁))
        #         @show "Broken at $t with $(norm(uₙ₊₁))"
        #         break
        #     end
        #     @show nei, te[nei], extrema(uₙ₊₁)
        # end
        # notify(uvis)
        @show "Done."
        io_enabled &&     vtk_save(pvd);

    end
    @time solve_direct(grid, topology, dh, refinement_cache, sync, 0.0001, false)

end
