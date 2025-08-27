using Ferrite
using Tensors
using Test
using Logging
using ForwardDiff
import SHA
using Random
using LinearAlgebra: norm2, svdvals, pinv, diag, mul!, muladd!
using SparseArrays
using StaticArrays
using OrderedCollections
using WriteVTK
import Metis
using HCubature: hcubature, hquadrature
using Unroll
using DataStructures
include("lts_utils.jl")
include("/home/amohamed/.julia/dev/Ferrite/src/Grid/Adaptivity/kopp.jl")
begin


    function analytical_solution(x, y, t; D=0.01, a=50)
        denominator = 1 + 4 * a * D * t
        return exp(-a * (x^2 + y^2) / denominator) / denominator
    end

    function solve_direct(grid::KoppGrid{Dim}, topology, dh, refinement_cache, sync, threshold, io_enabled::Bool=false) where Dim
        pvd = io_enabled ? paraview_collection("heat-benchmark-1-lts") : nothing
        T_end = 2.0
        ndofs_cell = 4
        u = sync.data_stores.solution_vector.data
        te = sync.data_stores.time_vector.data
        u_new = copy(u)
        u_prev = copy(u)
        error_vector = sync.data_stores.error_vector.data
        dofs_map = sync.data_stores.dofs_map
        dofs_temp_storage = zeros(Int, ndofs_cell)
        error_for_analysis = Tuple{Float64,Float64,Float64}[]
        for T in 0:0.1:T_end
            dt = Dict{Int,Float64}()
            nMinvK = Dict{Int,Matrix{Float64}}()
            te_pq = PriorityQueue{Int,Float64}()
            @time "prestep" @views for cc in CellIterator(grid)
                grid.kopp_cells[cellid(cc)].isleaf || continue
                offset_a = dofs_map.offsets[cellid(cc)]
                end_idx = cellid(cc) == length(dofs_map.offsets) ? length(dofs_map.dofs) : dofs_map.offsets[cellid(cc)+1] - 1
                K = sync.data_stores.assembled_stiffness_matrix.data[1:ndofs_cell, offset_a:end_idx]
                # display(sync.data_stores.cell_mass_matrix.data[:, :, cellid(cc)])
                Minv = inv((sync.data_stores.cell_mass_matrix.data[:, :, cellid(cc)]))
                celldofs!(dofs_temp_storage, dh, cellid(cc))
                nMinvK[cellid(cc)] = Minv * -K
                eig = maximum(sum(nMinvK[cellid(cc)], dims=2) + abs.(diag(nMinvK[cellid(cc)])) - diag(nMinvK[cellid(cc)]))
                # dt[cellid(cc)] = (1 / maximum(svdvals(nMinvK[cellid(cc)])))
                dt[cellid(cc)] = (1 / eig)
                enqueue!(te_pq, cellid(cc), te[cellid(cc)])
            end
            teprev = copy(te)
            u_prev = copy(u)
            if io_enabled
                dt_vec = zeros(length(te))
                @views for cc in CellIterator(grid)
                    grid.kopp_cells[cellid(cc)].isleaf || continue
                    offset_a = dofs_map.offsets[cellid(cc)]
                    end_idx = cellid(cc) == length(dofs_map.offsets) ? length(dofs_map.dofs) : dofs_map.offsets[cellid(cc)+1] - 1
                    K = sync.data_stores.assembled_stiffness_matrix.data[1:ndofs_cell, offset_a:end_idx]
                    Minv = inv((sync.data_stores.cell_mass_matrix.data[:, :, cellid(cc)]))
                    celldofs!(dofs_temp_storage, dh, cellid(cc))
                    nMinvK[cellid(cc)] = Minv * -K
                    dt_vec[cellid(cc)] = (1 / maximum(svdvals(nMinvK[cellid(cc)])))
                end
                fgrid, leaves = to_ferrite_grid(grid, true)
                dh2 = deepcopy(sync.dh)
                resize!(dh2.grid.cells, length(fgrid.cells))
                empty!(dh2.subdofhandlers[1].cellset)
                push!(dh2.subdofhandlers[1].cellset, (1:length(fgrid.cells))...)
                dh2.grid.cells .= fgrid.cells
                resize!(dh2.grid.nodes, length(fgrid.nodes))
                dh2.grid.nodes .= fgrid.nodes
                VTKGridFile("heat-benchmark-1-lts/heat-$T.vtu", fgrid; write_discontinuous=true) do vtk
                    write_solution(vtk, dh2, u, "_")
                    write_cell_data(vtk, error_vector[leaves], "__")
                    write_cell_data(vtk, dt_vec[leaves], "dt")
                    write_cell_data(vtk, [x.isleaf ? 1 : 0 for x in grid.kopp_cells], "isleaf")
                    pvd[T] = vtk
                end

            end
            # @views while any(t_el -> t_el < T, te)
            uhelp = Float64[]
            @time "iterating" @views @inbounds while last(peek(te_pq)) < T
                # u_new .= 0.0

                # error_vector .= 0.0

                nei = dequeue!(te_pq)

                if !grid.kopp_cells[nei].isleaf
                    te[nei] = Inf
                    continue
                end
                t = te[nei]
                Δt = dt[nei]
                if t + Δt > T
                    Δt = T - t
                end

                # @views for cc in CellIterator(grid)
                grid.kopp_cells[nei].isleaf || continue
                offset_a = dofs_map.offsets[nei]
                end_idx = nei == length(dofs_map.offsets) ? length(dofs_map.dofs) : dofs_map.offsets[nei+1] - 1
                dofs_a = dofs_map.dofs[offset_a:end_idx]
                nMinvK_local = nMinvK[nei]
                celldofs!(dofs_temp_storage, dh, nei)
                # @show Δt (1 / maximum(svdvals(Minv * (-K))))
                resize!(uhelp, size(nMinvK_local, 2))
                uhelp[1:ndofs_cell] .= u[dofs_a[1:ndofs_cell]]
                neighbor_idx = 1
                for lfi ∈ 1:2*Dim
                    for face_neighbor ∈ getneighborhood(topology, FacetIndex(nei, lfi))
                        dofs_neighbor_range = ndofs_cell+1+(neighbor_idx-1)*ndofs_cell:ndofs_cell+neighbor_idx*ndofs_cell
                        uhelp[dofs_neighbor_range] .= interpolate_linear.(teprev[face_neighbor[1]], te[face_neighbor[1]], u_prev[dofs_a[dofs_neighbor_range]], u[dofs_a[dofs_neighbor_range]], t)
                        # @show uhelp[dofs_neighbor_range] teprev[face_neighbor[1]] te[face_neighbor[1]] u_prev[dofs_a[dofs_neighbor_range]] u[dofs_a[dofs_neighbor_range]] t
                        # norm(uhelp[dofs_neighbor_range])>2.5 && return
                        neighbor_idx += 1
                    end
                end
                # if (T == 0.2) && nei == 61
                #     @show dofs_a getneighborhood(topology, FacetIndex(nei, 2))
                #     @show dofs_a getneighborhood(topology, FacetIndex(61, 4))
                #     return
                # end
                u_prev[dofs_a[1:ndofs_cell]] .= u[dofs_a[1:ndofs_cell]]

                mul!(u[dofs_temp_storage], nMinvK_local, uhelp, Δt, 1.0)
                # end
                # u[dofs_temp_storage] .= u_new[dofs_temp_storage]

                teprev[nei] = te[nei]
                te[nei] += Δt
                enqueue!(te_pq, nei, te[nei])

            end

            @time "error" update!(sync.data_stores.error_vector, sync.data_stores_prev.error_vector, sync.data_stores.solution_vector, grid, topology, dh, refinement_cache, sync.values_cache, 4)

            Base.resize!(sync.data_stores_prev.error_vector.data, size(sync.data_stores.error_vector.data)...)
            sync.data_stores_prev.error_vector.data .= sync.data_stores.error_vector.data
            if norm(u) > 1e3
                vtk_save(pvd)
                @show "Broken at $T"
                fgrid = to_ferrite_grid(grid)
                dh2 = deepcopy(sync.dh)
                resize!(dh2.grid.cells, length(fgrid.cells))
                dh2.grid.cells .= fgrid.cells
                resize!(dh2.grid.nodes, length(fgrid.nodes))
                dh2.grid.nodes .= fgrid.nodes
                VTKGridFile("amr-$T-broken.vtu", fgrid; write_discontinuous=true) do vtk
                    write_solution(vtk, dh2, u, "_")
                    write_cell_data(vtk, error_vector, "__")
                end
                break
            end
            errors = check_and_compute_convergence_norms(grid, dh, u, sync.values_cache.cell_values, T)
            push!(error_for_analysis, errors)
            needs_refinement = true
            refinement_iteration = 0
            @time "amr" begin
                # while needs_refinement == true
                needs_refinement = false
                refinement_set = OrderedSet{CellIndex}()
                coarsening_set = OrderedSet{CellIndex}()
                sizehint!(refinement_set, length(grid.kopp_cells))
                sizehint!(coarsening_set, length(grid.kopp_cells))
                coarsening_vector = zeros(Int, length(grid.kopp_cells))
                @inbounds for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                    grid.kopp_cells[cellid(cc)].isleaf || continue
                    if sqrt(error_vector[cellid(cc)]) > threshold
                        get_refinement_level(grid.kopp_cells[cellid(cc)]) >= 5 && continue
                        push!(refinement_set, CellIndex(cellid(cc)))
                        needs_refinement = true
                    end
                end
                isempty(refinement_set) || refine!(grid, topology, refinement_cache, sync, refinement_set)
                coarsening_vector = zeros(Int, length(grid.kopp_cells))
                for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                    grid.kopp_cells[cellid(cc)].isleaf || continue
                    if threshold / 10 > sqrt(error_vector[cellid(cc)])
                        parent = grid.kopp_cells[cellid(cc)].parent
                        parent <= 0 && continue
                        coarsening_vector[parent] += 1
                    end
                end
                for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                    grid.kopp_cells[cellid(cc)].isleaf && continue
                    coarsening_vector[cellid(cc)] < 4 && continue
                    push!(coarsening_set, CellIndex(cellid(cc)))
                    needs_refinement = true
                end
                isempty(coarsening_set) || coarsen!(grid, topology, refinement_cache, sync, coarsening_set)
                refinement_iteration += 1
                # refinement_iteration == 2 && break
                @warn T
                # end
            end
            resize!(u_new, length(u))

        end
        # notify(uvis)
        io_enabled && vtk_save(pvd)
        return error_for_analysis
    end

    function time_once_threshod(threshold)
        refshape = RefQuadrilateral
        grid = generate_grid(KoppCell{2,Int}, (10, 10), Ferrite.Vec((-1.0, -1.0)), Ferrite.Vec((1.0, 1.0)))

        topology = KoppTopology(grid)

        ip = DiscontinuousLagrange{refshape,1}()
        qr = QuadratureRule{refshape}(2)

        face_qr = FacetQuadratureRule{refshape}(2)
        cellvalues = CellValues(qr, ip)
        facetvalues = FacetValues(face_qr, ip)
        interfacevalues = InterfaceValues(face_qr, ip)

        dh = DofHandler(grid.base_grid)
        add!(dh, :u, ip)
        close!(dh)

        lts_values = ValuesCache(
            cellvalues,
            facetvalues,
            interfacevalues
        )

        refinement_cache = KoppRefinementCache(grid, topology)
        sync = LTSAMRSynchronizer(grid, dh, lts_values, refinement_cache, topology, 0.1)
        u = sync.data_stores.solution_vector.data
        Ferrite.apply_analytical!(u, dh, :u, x -> analytical_solution(x[1], x[2], 0.0))
        needs_refinement = true
        update!(sync.data_stores.error_vector, sync.data_stores_prev.error_vector, sync.data_stores.solution_vector, grid, topology, dh, refinement_cache, sync.values_cache, 4)
        error_vector = sync.data_stores.error_vector.data
        while needs_refinement == true
            needs_refinement = false
            fgrid, ids = to_ferrite_grid(grid)
            dh2 = deepcopy(sync.dh)
            resize!(dh2.grid.cells, length(fgrid.cells))
            dh2.grid.cells .= fgrid.cells
            resize!(dh2.grid.nodes, length(fgrid.nodes))
            dh2.grid.nodes .= fgrid.nodes
            Ferrite.apply_analytical!(u, dh2, :u, x -> analytical_solution(x[1], x[2], 0.0))
            needs_refinement = false
            refinement_set = OrderedSet{CellIndex}()
            coarsening_set = OrderedSet{CellIndex}()
            sizehint!(refinement_set, length(grid.kopp_cells))
            sizehint!(coarsening_set, length(grid.kopp_cells))
            coarsening_vector = zeros(Int, length(grid.kopp_cells))
            @inbounds for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                grid.kopp_cells[cellid(cc)].isleaf || continue
                if sqrt(error_vector[cellid(cc)]) > threshold
                    get_refinement_level(grid.kopp_cells[cellid(cc)]) >= 5 && continue
                    push!(refinement_set, CellIndex(cellid(cc)))
                    needs_refinement = true
                end
            end
            isempty(refinement_set) || refine!(grid, topology, refinement_cache, sync, refinement_set)
            coarsening_vector = zeros(Int, length(grid.kopp_cells))
            for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                grid.kopp_cells[cellid(cc)].isleaf || continue
                if threshold / 10 > sqrt(error_vector[cellid(cc)])
                    parent = grid.kopp_cells[cellid(cc)].parent
                    parent <= 0 && continue
                    coarsening_vector[parent] += 1
                end
            end
            for cc in CellIterator(grid, 1:length(grid.kopp_cells))
                grid.kopp_cells[cellid(cc)].isleaf && continue
                coarsening_vector[cellid(cc)] < 4 && continue
                push!(coarsening_set, CellIndex(cellid(cc)))
                needs_refinement = true
            end
            isempty(coarsening_set) || coarsen!(grid, topology, refinement_cache, sync, coarsening_set)
        end


        @time error_for_analysis = solve_direct(grid, topology, dh, refinement_cache, sync, threshold, false)
        return error_for_analysis
    end
    using GLMakie
    errors_viz = Float64[]
    thresholds = 0.01:0.001:0.1
    for threshold in thresholds
        t_error = 0.0
        error_for_analysis = time_once_threshod(threshold)
        dt = 0.1
        for (i, t) in enumerate(0.0:dt:2.0)
            t_error += error_for_analysis[i][1] * dt
        end
        push!(errors_viz, t_error)
    end
end
