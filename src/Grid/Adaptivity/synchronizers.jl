abstract type AbstractAMRSynchronizer end

struct NullAMRSynchronizer <: AbstractAMRSynchronizer end

struct LTSAMRSynchronizer{V,VCT<:ValuesCache,DHT<:DofHandler} <: AbstractAMRSynchronizer
    # Values
    values_cache::VCT
    # DofHandler
    dh::DHT
    celldofs_prev::Vector{Int}
    cell_dofs_offset_prev::Vector{Int}
    cell_to_subdofhandler_prev::Vector{Int}
    # Mass matrices
    M_cell_matrices::ElasticArray{Float64,3,2,V}
    # Stiffness matrices
    K_cell_matrices::ElasticArray{Float64,3,2,V}
    K_interface_matrices::ElasticArray{Float64,3,2,V}
    # Source vector
    N_cell_vectors::Vector{Float64}
    # Solution vector
    u::Vector{Float64}
    # Indexing
    interface_matrix_index::Vector{Int}
end

function sync_amr_refinement_forward!(grid::KoppGrid, sync::LTSAMRSynchronizer, refinement_cache::KoppRefinementCache, n_refined_cells::Int, n_neighborhoods::Int)
    Dim = Ferrite.getspatialdim(grid.base_grid)
    new_length = last(size(sync.M_cell_matrices))+(2^Dim)*n_refined_cells
    resize!(sync.M_cell_matrices, (size(sync.M_cell_matrices)[1], size(sync.M_cell_matrices)[2], new_length))
    resize!(sync.K_cell_matrices, (size(sync.K_cell_matrices)[1], size(sync.K_cell_matrices)[2], new_length))
    # resize!(sync.K_interface_matrix_index, n_neighborhoods)
    resize!(sync.K_interface_matrices, (size(sync.K_interface_matrices)[1], size(sync.K_interface_matrices)[2], n_neighborhoods ÷ 2))
    # resize!(sync.ansatz_isactive, new_length * sync.dh.subdofhandlers[1].ndofs_per_cell)
    resize!(sync.dh.cell_dofs, new_length * sync.dh.subdofhandlers[1].ndofs_per_cell)
    resize!(sync.dh.cell_dofs_offset, new_length)
    resize!(sync.dh.cell_to_subdofhandler, new_length)

    update_dofs!(grid, refinement_cache, sync.dh, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev)

    sync.interface_matrix_index .= 0
end

function sync_amr_refinement_backward!(sync::LTSAMRSynchronizer)
    copy!(sync.celldofs_prev, sync.dh.cell_dofs)
    copy!(sync.cell_dofs_offset_prev, sync.dh.cell_dofs_offset)
    copy!(sync.cell_to_subdofhandler_prev, sync.dh.cell_to_subdofhandler)

end

function LTSAMRSynchronizer(grid::KoppGrid, dh::Ferrite.AbstractDofHandler, kopp_values::ValuesCache, refinement_cache::KoppRefinementCache, topology::KoppTopology{NFacets}, u = zeros(ndofs(dh))) where NFacets
    @assert Ferrite.isclosed(dh) "DofHandler must be closed"
    ndofs_cell = ndofs_per_cell(dh)
    ref_shape = getrefshape(grid.kopp_cells[1])
    ip = DiscontinuousLagrange{ref_shape, 1}()
    qr = QuadratureRule{ref_shape}(1)
    geo_mapping = Ferrite.GeometryMapping{1}(Float64, ip, qr)
    # M
    M_cell_matrices = ElasticArray(zeros(Float64, ndofs_cell, ndofs_cell, length(refinement_cache.children_updated_indices)))
    for cell_cache in Ferrite.CellIterator(grid)
        reinit!(kopp_values.cell_values, cell_cache)
        cell_idx = cell_cache.cellid
        assemble_element_matrix!((@view M_cell_matrices[:,:,cell_idx]), kopp_values)
    end
    # K_cell_matrices
    K_cell_matrices = ElasticArray(zeros(Float64, ndofs_cell, ndofs_cell, length(refinement_cache.children_updated_indices)))
    for cell_cache in Ferrite.CellIterator(grid)
        reinit!(kopp_values.cell_values, cell_cache)
        cell_idx = cell_cache.cellid
        assemble_element_matrix!((@view K_cell_matrices[:,:,cell_idx]), kopp_values)
    end
    K_interface_matrices = ElasticArray(zeros(Float64, 2 * ndofs_cell, 2 * ndofs_cell, length(topology.neighbors)÷2))
    # Construct interface matrix indices
    interface_matrix_indices = Vector{Int}(undef, count(topology.cell_facet_neighbors_length .!= 0 ))
    interface_index = 1
    ii = Ferrite.InterfaceIterator(grid, topology)
    for (facet_idx, nneighbors, ninterfaces) in ii
        interface_cache = ii.cache
        reinit!(kopp_values.interface_values, interface_cache, topology)
        assemble_interface_matrix!((@view K_interface_matrices[:,:,interface_index]), kopp_values)
        interface_index += 1
    end
    display(K_interface_matrices)

    N_cell_vectors = zero(u)
    new_offset = 1
    for (idx, offset) in pairs(IndexCartesian(), topology.cell_facet_neighbors_offset)
        offset == 0 && continue
        neighbor_cell = topology.neighbors[offset][1]
        neighbor_facet = topology.neighbors[offset][2]
        current_cell = idx[2]
        if(neighbor_cell < current_cell)
            interface_matrix_indices[offset] = interface_matrix_indices[topology.cell_facet_neighbors_offset[neighbor_facet, neighbor_cell]]
            continue
        end
        interface_matrix_indices[offset] = new_offset
        new_offset += 1
    end
    ansatz_isactive = trues(ndofs(dh))
    return LTSAMRSynchronizer(kopp_values,
                    dh,
                    copy(dh.cell_dofs),
                    copy(dh.cell_dofs_offset),
                    copy(dh.cell_to_subdofhandler),
                    M_cell_matrices,
                    K_cell_matrices,
                    K_interface_matrices,
                    N_cell_vectors,
                    u,
                    interface_matrix_indices)
end

struct KoppSynchronizer end
