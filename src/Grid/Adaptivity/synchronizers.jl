using Unrolled

abstract type AbstractAMRSynchronizer end

struct NullAMRSynchronizer <: AbstractAMRSynchronizer end

abstract type AMRDataContainer end

abstract type AMRCellData{C} <: AMRDataContainer end

struct CellMassMatrix{C} <: AMRCellData{C}
    data::C
end

struct CellStiffnessMatrix{C} <: AMRCellData{C}
    data::C
end

abstract type AMRInterfaceData{C} <: AMRDataContainer end

struct InterfaceStiffnessMatrix{C} <: AMRInterfaceData{C}
    data::C
end

abstract type AMRdofwiseVector{C} <: AMRDataContainer end

struct SolutionVector{C} <: AMRdofwiseVector{C}
    data::C
end

struct TimeVector{C} <: AMRCellData{C}
    data::C
end

function assemble_element_matrix!(K_cell_matrices, kopp_values, cell_idx) error("FUCK MICROSOFT") end

function Base.resize!(data_store::TimeVector, data_store_prev::TimeVector, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store.data, ncells)
    # data_store.data .= zero(eltype(data_store.data))
end

function Base.resize!(data_store::AMRdofwiseVector, data_store_prev::AMRdofwiseVector, ncells, ninterfaces, ndofs)
    Base.resize!(data_store_prev.data, size(data_store.data)...)
    data_store_prev.data .= data_store.data
    Base.resize!(data_store.data, ndofs)
    # data_store.data .= zero(eltype(data_store.data))
end

function Base.resize!(data_store::AMRCellData{<:ElasticArray}, data_store_prev::AMRCellData{<:ElasticArray}, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store.data, size(data_store.data)[1], size(data_store.data)[2], ncells)
    data_store.data .= zero(eltype(data_store.data))
end

function Base.resize!(data_store::AMRInterfaceData{<:ElasticArray}, data_store_prev::AMRInterfaceData{<:ElasticArray}, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store.data, size(data_store.data)[1], size(data_store.data)[2], ninterfaces)
    data_store.data .= zero(eltype(data_store.data))
end

function update!(data_store::TimeVector, data_store_prev::TimeVector, grid, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int})

end

function update!(data_store::AMRdofwiseVector, data_store_prev::AMRdofwiseVector, grid::KoppGrid{Dim}, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int}) where Dim
    data_store.data .= 0.0
    sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[1]]
    empty!(sdh.cellset)
    push!(sdh.cellset, (1:length(grid.kopp_cells))...)
    @views for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        new == 0 && continue
        new_offset = dh.cell_dofs_offset[new]
        new_sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[new]]
        # TODO: do only once
        new_dofs = dh.cell_dofs[new_offset : new_offset + new_sdh.ndofs_per_cell - 1]
        old_offset = cell_dofs_offset_prev[old]
        old_sdh = dh.subdofhandlers[cell_to_subdofhandler_prev[old]]
        old_dofs = celldofs_prev[old_offset : old_offset + old_sdh.ndofs_per_cell - 1]
        data_store.data[new_dofs] .= data_store_prev.data[old_dofs]
    end

    for cell_idx in 1:length(grid.kopp_cells)
        cell = getcells(grid, cell_idx)
        # Calculate new interfaces
        if refinement_cache.marked_for_refinement[cell_idx]
            parent_offset = dh.cell_dofs_offset[cell_idx]
            parent_sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx]]
            parent_dofs = @view dh.cell_dofs[parent_offset : parent_offset + parent_sdh.ndofs_per_cell - 1]
            for new_seq in 1 : 2^Dim
                child_idx = cell_idx + new_seq
                # DofHandler add new dofs and offsets
                # For DG only
                offset = dh.cell_dofs_offset[child_idx]
                sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[child_idx]]
                dofs = @view dh.cell_dofs[offset : offset + sdh.ndofs_per_cell - 1]
                # for dof in dofs
                #     @info new_seq dof data_store.data[dof]
                # end
                # interpolate_solution!(Ferrite.RefHypercube{Dim}, data_store.data, dofs, parent_dofs, new_seq)
            end
        end
        if refinement_cache.marked_for_coarsening[cell_idx]
        end
    end
end

function update!(data_store::AMRCellData{<:ElasticArray}, data_store_prev::AMRCellData{<:ElasticArray}, grid, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int})

end

function update!(data_store::AMRInterfaceData{<:ElasticArray}, data_store_prev::AMRInterfaceData{<:ElasticArray}, grid, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int})

end

struct LTSAMRSynchronizer{DataStores,VCT<:ValuesCache,DHT<:DofHandler} <: AbstractAMRSynchronizer
    # Values
    values_cache::VCT
    # DofHandler
    dh::DHT
    celldofs_prev::Vector{Int}
    cell_dofs_offset_prev::Vector{Int}
    cell_to_subdofhandler_prev::Vector{Int}
    # Where matrices and vectors are stored
    data_stores_prev::DataStores # Tuple of matrices and/or vectors wrapped in AMRDataContainer
    data_stores::DataStores # Tuple of matrices and/or vectors wrapped in AMRDataContainer
    # Indexing
    interface_matrix_index::Vector{Int}
    # scalars
    Δt::Float64
end

function sync_amr_refinement_forward!(grid::KoppGrid, sync::LTSAMRSynchronizer, refinement_cache::KoppRefinementCache, n_refined_cells::Int, n_neighborhoods::Int)
    Dim = Ferrite.getspatialdim(grid.base_grid)
    new_length = Ferrite.getncells(grid)+(2^Dim)*n_refined_cells

    unrolled_foreach((data_store_) -> Base.resize!(data_store_[1],  data_store_[2], new_length, n_neighborhoods ÷ 2, maximum(sync.dh.cell_dofs) + n_refined_cells * (2^Dim - 1) * sync.dh.subdofhandlers[1].ndofs_per_cell), zip(sync.data_stores, sync.data_stores_prev))

    resize!(sync.celldofs_prev, length(sync.dh.cell_dofs))
    resize!(sync.cell_dofs_offset_prev, length(sync.dh.cell_dofs_offset))
    resize!(sync.cell_to_subdofhandler_prev, length(sync.dh.cell_to_subdofhandler))

    copy!(sync.celldofs_prev, (sync.dh.cell_dofs))
    copy!(sync.cell_dofs_offset_prev, (sync.dh.cell_dofs_offset))
    copy!(sync.cell_to_subdofhandler_prev, (sync.dh.cell_to_subdofhandler))

    resize!(sync.dh.cell_dofs, new_length * sync.dh.subdofhandlers[1].ndofs_per_cell)
    resize!(sync.dh.cell_dofs_offset, new_length)
    resize!(sync.dh.cell_to_subdofhandler, new_length)

    update_dofs!(grid, refinement_cache, sync.dh, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev)

    sync.interface_matrix_index .= 0
end

function sync_amr_coarsening_forward!(grid::KoppGrid, sync::LTSAMRSynchronizer, refinement_cache::KoppRefinementCache, n_coarsened_cells::Int, n_neighborhoods::Int)
    Dim = Ferrite.getspatialdim(grid.base_grid)
    new_length = Ferrite.getncells(grid)-(2^Dim)*n_coarsened_cells

    unrolled_foreach((data_store_) -> Base.resize!(data_store_[1],  data_store_[2], new_length, n_neighborhoods ÷ 2, maximum(sync.dh.cell_dofs) - n_coarsened_cells * (2^Dim - 1) * sync.dh.subdofhandlers[1].ndofs_per_cell), zip(sync.data_stores, sync.data_stores_prev))

    resize!(sync.celldofs_prev, length(sync.dh.cell_dofs))
    resize!(sync.cell_dofs_offset_prev, length(sync.dh.cell_dofs_offset))
    resize!(sync.cell_to_subdofhandler_prev, length(sync.dh.cell_to_subdofhandler))

    copy!(sync.celldofs_prev, (sync.dh.cell_dofs))
    copy!(sync.cell_dofs_offset_prev, (sync.dh.cell_dofs_offset))
    copy!(sync.cell_to_subdofhandler_prev, (sync.dh.cell_to_subdofhandler))

    resize!(sync.dh.cell_dofs, new_length * sync.dh.subdofhandlers[1].ndofs_per_cell)
    resize!(sync.dh.cell_dofs_offset, new_length)
    resize!(sync.dh.cell_to_subdofhandler, new_length)

    update_dofs!(grid, refinement_cache, sync.dh, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev)

    sync.interface_matrix_index .= 0
end

function sync_amr_refinement_backward!(sync::LTSAMRSynchronizer, refinement_cache::KoppRefinementCache, grid, topology)

    # for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
    #     # copy old cells
    #     kopp_cache.cell_matrices[:,:,new] .= @view kopp_cache.cell_matrices[:,:,old]
    #     old != new && (kopp_cache.cell_matrices[:,:,new] .= 0.0)
    # end
    # for cell_cache in Ferrite.CellIterator(grid, UpdateFlags(false, true, true))
    #     cell_idx = cell_cache.cellid
    #     cell = getcells(grid, cell_idx)
    #     # Calculate new interfaces
    #     if refinement_cache.marked_for_refinement[cell_idx]
    #         parent_offset = dh.cell_dofs_offset[cell_idx]
    #         parent_sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx]]
    #         parent_dofs = @view dh.cell_dofs[parent_offset : parent_offset + parent_sdh.ndofs_per_cell - 1]
    #         for new_seq in 1 : 2^Dim
    #             child_idx = cell_idx + new_seq
    #             # DofHandler add new dofs and offsets
    #             # For DG only
    #             offset = dh.cell_dofs_offset[child_idx]
    #             sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[child_idx]]
    #             dofs = @view dh.cell_dofs[offset : offset + sdh.ndofs_per_cell - 1]
    #             interpolate_solution!(Ferrite.RefHypercube{Dim}, kopp_cache.u, dofs, parent_dofs, new_seq)
    #         end
    #     end
    #     if all(x -> x == 0.0, @view kopp_cache.cell_matrices[:,:,cell_idx])
    #         reinit!(kopp_values.cell_values, cell_cache)
    #         assemble_element_matrix!((@view kopp_cache.cell_matrices[:,:,cell_idx]), kopp_values)
    #     end
    # end
    # ninterfaces = 1
    # facet_idx = 1
    # nneighbors = 1
    # interface_index = 1
    # ii = Ferrite.InterfaceIterator(dh, topology)
    # for (facet_idx, nneighbors, ninterfaces) in ii
    #     interface_cache = ii.cache
    #     _facet_idx = nneighbors == 1 ? facet_idx - 1 : facet_idx
    #     cell_idx = (_facet_idx - 1) ÷ (2*Dim) + 1
    #     facet_a = (_facet_idx - 1) % (2*Dim) + 1
    #     neighbor = FacetIndex(interface_cache.b.cc.cellid, interface_cache.b.current_facet_id)
    #     kopp_cache.interface_matrix_index[topology.cell_facet_neighbors_offset[facet_a, cell_idx] + nneighbors - 1] = ninterfaces
    #     rng  = topology.cell_facet_neighbors_offset[neighbor[2], neighbor[1]] : topology.cell_facet_neighbors_offset[neighbor[2], neighbor[1]] + topology.cell_facet_neighbors_length[neighbor[2], neighbor[1]] - 1
    #     # @assert any(x -> x == 0, @view kopp_cache.interface_matrix_index[rng])
    #     for j in rng
    #         kopp_cache.interface_matrix_index[j] != 0 && continue
    #         kopp_cache.interface_matrix_index[j] = ninterfaces
    #         break
    #     end
    #     reinit!(kopp_values.interface_values, interface_cache, topology)
    #     assemble_interface_matrix!((@view kopp_cache.interface_matrices[:,:,interface_index]), kopp_values)
    #     interface_index += 1
    # end
    unrolled_foreach((data_store_) -> update!(data_store_[1],  data_store_[2], grid, topology, sync.dh, refinement_cache, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev), zip(sync.data_stores, sync.data_stores_prev))

    copy!(sync.celldofs_prev, sync.dh.cell_dofs)
    copy!(sync.cell_dofs_offset_prev, sync.dh.cell_dofs_offset)
    copy!(sync.cell_to_subdofhandler_prev, sync.dh.cell_to_subdofhandler)

end

function sync_amr_coarsening_backward!(sync::LTSAMRSynchronizer, refinement_cache::KoppRefinementCache, grid, topology)

    # for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
    #     # copy old cells
    #     kopp_cache.cell_matrices[:,:,new] .= @view kopp_cache.cell_matrices[:,:,old]
    #     old != new && (kopp_cache.cell_matrices[:,:,new] .= 0.0)
    # end
    # for cell_cache in Ferrite.CellIterator(grid, UpdateFlags(false, true, true))
    #     cell_idx = cell_cache.cellid
    #     cell = getcells(grid, cell_idx)
    #     # Calculate new interfaces
    #     if refinement_cache.marked_for_refinement[cell_idx]
    #         parent_offset = dh.cell_dofs_offset[cell_idx]
    #         parent_sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx]]
    #         parent_dofs = @view dh.cell_dofs[parent_offset : parent_offset + parent_sdh.ndofs_per_cell - 1]
    #         for new_seq in 1 : 2^Dim
    #             child_idx = cell_idx + new_seq
    #             # DofHandler add new dofs and offsets
    #             # For DG only
    #             offset = dh.cell_dofs_offset[child_idx]
    #             sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[child_idx]]
    #             dofs = @view dh.cell_dofs[offset : offset + sdh.ndofs_per_cell - 1]
    #             interpolate_solution!(Ferrite.RefHypercube{Dim}, kopp_cache.u, dofs, parent_dofs, new_seq)
    #         end
    #     end
    #     if all(x -> x == 0.0, @view kopp_cache.cell_matrices[:,:,cell_idx])
    #         reinit!(kopp_values.cell_values, cell_cache)
    #         assemble_element_matrix!((@view kopp_cache.cell_matrices[:,:,cell_idx]), kopp_values)
    #     end
    # end
    # ninterfaces = 1
    # facet_idx = 1
    # nneighbors = 1
    # interface_index = 1
    # ii = Ferrite.InterfaceIterator(dh, topology)
    # for (facet_idx, nneighbors, ninterfaces) in ii
    #     interface_cache = ii.cache
    #     _facet_idx = nneighbors == 1 ? facet_idx - 1 : facet_idx
    #     cell_idx = (_facet_idx - 1) ÷ (2*Dim) + 1
    #     facet_a = (_facet_idx - 1) % (2*Dim) + 1
    #     neighbor = FacetIndex(interface_cache.b.cc.cellid, interface_cache.b.current_facet_id)
    #     kopp_cache.interface_matrix_index[topology.cell_facet_neighbors_offset[facet_a, cell_idx] + nneighbors - 1] = ninterfaces
    #     rng  = topology.cell_facet_neighbors_offset[neighbor[2], neighbor[1]] : topology.cell_facet_neighbors_offset[neighbor[2], neighbor[1]] + topology.cell_facet_neighbors_length[neighbor[2], neighbor[1]] - 1
    #     # @assert any(x -> x == 0, @view kopp_cache.interface_matrix_index[rng])
    #     for j in rng
    #         kopp_cache.interface_matrix_index[j] != 0 && continue
    #         kopp_cache.interface_matrix_index[j] = ninterfaces
    #         break
    #     end
    #     reinit!(kopp_values.interface_values, interface_cache, topology)
    #     assemble_interface_matrix!((@view kopp_cache.interface_matrices[:,:,interface_index]), kopp_values)
    #     interface_index += 1
    # end
    unrolled_foreach((data_store_) -> update!(data_store_[1],  data_store_[2], grid, topology, sync.dh, refinement_cache, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev), zip(sync.data_stores, sync.data_stores_prev))

    copy!(sync.celldofs_prev, sync.dh.cell_dofs)
    copy!(sync.cell_dofs_offset_prev, sync.dh.cell_dofs_offset)
    copy!(sync.cell_to_subdofhandler_prev, sync.dh.cell_to_subdofhandler)

end

function assemble_element_matrix!(K::CellMassMatrix, kopp_values::ValuesCache, cell_idx)
    Ke = @view K.data[:,:,cell_idx]
    cv = kopp_values.cell_values
    n_basefuncs = getnbasefunctions(cv)
    for q_point in 1:getnquadpoints(cv.qr)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs
            δu  = Ferrite.shape_value(cv, q_point, i)
            for j in 1:n_basefuncs
                u = Ferrite.shape_value(cv, q_point, i)
                Ke[i, j] += (δu ⋅ u) * dΩ
            end
        end
    end
    return nothing
end
function assemble_element_matrix!(K::CellStiffnessMatrix, kopp_values::ValuesCache, cell_idx)
    Ke = @view K.data[:,:,cell_idx]
    cv = kopp_values.cell_values
    n_basefuncs = getnbasefunctions(cv)
    for q_point in 1:getnquadpoints(cv)
        ## Get the quadrature weight
        dΩ = getdetJdV(cv, q_point)
        ## Loop over test shape functions
        for i in 1:n_basefuncs
            ∇δu = shape_gradient(cv, q_point, i)
            ## Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cv, q_point, j)
                ## Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return nothing
end
function assemble_element_matrix!(K::InterfaceStiffnessMatrix, kopp_values::ValuesCache, interface_index, μ::Float64 = 10.)
    Ki = @view K.data[:,:,interface_index]
    iv = kopp_values.interface_values
    for q_point in 1:getnquadpoints(iv)
        # Get the normal to facet A
        normal = getnormal(iv, q_point)
        # Get the quadrature weight
        dΓ = getdetJdV(iv, q_point)
        # Loop over test shape functions
        for i in 1:getnbasefunctions(iv)
            # Multiply the jump by the negative normal to get the definition from the theory section.
            δu_jump = shape_value_jump(iv, q_point, i) * (-normal)
            ∇δu_avg = shape_gradient_average(iv, q_point, i)
            # Loop over trial shape functions
            for j in 1:getnbasefunctions(iv)
                # Multiply the jump by the negative normal to get the definition from the theory section.
                u_jump = shape_value_jump(iv, q_point, j) * (-normal)
                ∇u_avg = shape_gradient_average(iv, q_point, j)
                # Add contribution to Ki
                Ki[i, j] += -(δu_jump ⋅ ∇u_avg + ∇δu_avg ⋅ u_jump)*dΓ +  μ * (δu_jump ⋅ u_jump) * dΓ
            end
        end
    end
end


function LTSAMRSynchronizer(grid::KoppGrid, dh::Ferrite.AbstractDofHandler, kopp_values::ValuesCache, refinement_cache::KoppRefinementCache, topology::KoppTopology{NFacets}, Δt = 0.1, _u = zeros(ndofs(dh))) where NFacets
    @assert Ferrite.isclosed(dh) "DofHandler must be closed"
    ndofs_cell = ndofs_per_cell(dh)
    # M
    M_cell_matrices = CellMassMatrix(ElasticArray(zeros(Float64, ndofs_cell, ndofs_cell, length(refinement_cache.old_cell_to_new_cell_map))))
    for cell_cache in Ferrite.CellIterator(grid)
        reinit!(kopp_values.cell_values, cell_cache)
        cell_idx = cell_cache.cellid
        assemble_element_matrix!(M_cell_matrices, kopp_values, cell_idx)
    end
    # K_cell_matrices
    K_cell_matrices = CellStiffnessMatrix(ElasticArray(zeros(Float64, ndofs_cell, ndofs_cell, length(refinement_cache.old_cell_to_new_cell_map))))
    for cell_cache in Ferrite.CellIterator(grid)
        reinit!(kopp_values.cell_values, cell_cache)
        cell_idx = cell_cache.cellid
        assemble_element_matrix!(K_cell_matrices, kopp_values, cell_idx)
    end
    K_interface_matrices = InterfaceStiffnessMatrix(ElasticArray(zeros(Float64, 2 * ndofs_cell, 2 * ndofs_cell, length(topology.neighbors)÷2)))
    # Construct interface matrix indices
    interface_matrix_indices = Vector{Int}(undef, count(topology.cell_facet_neighbors_length .!= 0 ))
    interface_index = 1
    ii = Ferrite.InterfaceIterator(grid, topology)
    for interface_cache in ii
        reinit!(kopp_values.interface_values, interface_cache, topology)
        assemble_element_matrix!(K_interface_matrices, kopp_values, interface_index)
        interface_index += 1
    end

    # N_cell_vectors = zero(u)
    # new_offset = 1
    # for (idx, offset) in pairs(IndexCartesian(), topology.cell_facet_neighbors_offset)
    #     offset == 0 && continue
    #     neighbor_cell = topology.neighbors[offset][1]
    #     neighbor_facet = topology.neighbors[offset][2]
    #     current_cell = idx[2]
    #     if(neighbor_cell < current_cell)
    #         interface_matrix_indices[offset] = interface_matrix_indices[topology.cell_facet_neighbors_offset[neighbor_facet, neighbor_cell]]
    #         continue
    #     end
    #     interface_matrix_indices[offset] = new_offset
    #     new_offset += 1
    # end
    u = SolutionVector(_u)
    t = TimeVector(zeros(getncells(grid)))
    ansatz_isactive = trues(ndofs(dh))
    return LTSAMRSynchronizer(kopp_values,
                    dh,
                    copy(dh.cell_dofs),
                    copy(dh.cell_dofs_offset),
                    copy(dh.cell_to_subdofhandler),
                    (M_cell_matrices, K_cell_matrices, K_interface_matrices, u, t),
                    (deepcopy(M_cell_matrices), deepcopy(K_cell_matrices), deepcopy(K_interface_matrices), deepcopy(u), deepcopy(t)),
                    interface_matrix_indices,
                    Δt)
end

struct KoppSynchronizer end
