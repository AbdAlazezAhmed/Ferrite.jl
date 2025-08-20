function __resize_marked_for_refinement!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache,
    cellset::AbstractSet{CellIndex}) where {Dim}
    n_refined_cells = length(cellset)
    new_length = length(grid.kopp_cells) - (2^Dim) * n_refined_cells
    resize!(refinement_cache.marked_for_coarsening, new_length)
    resize!(refinement_cache.new_cell_to_old_cell_map, new_length)

    refinement_cache.marked_for_coarsening[length(grid.kopp_cells) + 1 : end] .= false
    refinement_cache.new_cell_to_old_cell_map .= 0
    return nothing
end

function __update_refinement_cache_isactive!(
    grid::KoppGrid{Dim},
    topology::KoppTopology, # Old topology
    refinement_cache::KoppRefinementCache,
    cellset::AbstractSet{CellIndex}
    ) where {Dim}
    n_refined_interfaces = 0
    for (i, cell) in enumerate(grid.kopp_cells)
        if CellIndex(i) ∈ cellset
            refinement_cache.old_cell_to_new_cell_map[i+(2^Dim)+1:end] .= (@view refinement_cache.old_cell_to_new_cell_map[i+(2^Dim)+1:end]) .- 2^Dim
            refinement_cache.old_cell_to_new_cell_map[i+1:i+(2^Dim)] .= 0
            refinement_cache.marked_for_coarsening[refinement_cache.old_cell_to_new_cell_map[i]] = true
            # kopp_cache.ansatz_isactive[@view dh.cell_dofs[dh.cell_dofs_offset[i]:dh.cell_dofs_offset[i]+dh.subdofhandlers[1].ndofs_per_cell-1]] .= true
            # kopp_cache.ansatz_isactive[dh.cell_dofs_offset[i]+dh.subdofhandlers[1].ndofs_per_cell : dh.cell_dofs_offset[i]+dh.subdofhandlers[1].ndofs_per_cell*(2^Dim) - 1] .= true
            for facet in 1:2Dim
                interface_offset = topology.cell_facet_neighbors_offset[facet, i]
                interface_length = topology.cell_facet_neighbors_length[facet, i]
                for interface_index_idx in interface_offset : interface_offset + interface_length - 1
                    interface_offset == 0 && continue
                    neighbor = topology.neighbors[interface_index_idx]
                    # Refine interface when the current cell is lower index than the neighbor only to avoid race conditions
                    get_refinement_level(grid.kopp_cells[neighbor.idx[1]]) == get_refinement_level(cell) + 1 && i > neighbor.idx[1] && continue
                    # Refine the interface only if the neighbor is not finer than the current cell
                    get_refinement_level(grid.kopp_cells[neighbor.idx[1]]) > get_refinement_level(cell) + 1 && continue
                    # TODO: Invalidate interface matrices
                    # refinement_cache.interfaces_updated_indices[interface_index_idx+1:end] .= (@view refinement_cache.interfaces_updated_indices[interface_index_idx+1:end]) .- (2^(Dim-1) - 1)
                    n_refined_interfaces += 1
                end
            end
        end
        new_idx = refinement_cache.old_cell_to_new_cell_map[i]
        iszero(new_idx) || (refinement_cache.new_cell_to_old_cell_map[new_idx] = i)
    end
    # old_length = length(refinement_cache.interfaces_data_updated_indices)
    # resize!(refinement_cache.interfaces_data_updated_indices, old_length - n_refined_interfaces * ((2^(Dim-1) - 1)) - 2Dim * length(cellset))
    # refinement_cache.interfaces_data_updated_indices .= 0
    # old_length = length(refinement_cache.interfaces_updated_indices)
    # resize!(refinement_cache.interfaces_updated_indices, old_length - n_refined_interfaces * ((2^(Dim-1) - 1)) - 2Dim * length(cellset))
    # refinement_cache.interfaces_data_updated_indices .= 0
    return nothing
end

function _count_neighbors_update_indexing!(
    grid::KoppGrid{Dim},
    topology::KoppTopology,
    refinement_cache::KoppRefinementCache) where {Dim}
    NFacets = 2 * Dim
    temp_cell_facet_neighbors_offset = topology.cell_facet_neighbors_offset_prev
    temp_cell_facet_neighbors_length = topology.cell_facet_neighbors_length_prev
    n_neighborhoods = 0
    ref_cell_neighborhood = get_ref_cell_neighborhood(getrefshape(getcelltype(grid.base_grid)))
    set = Set{FacetIndex}()
    sizehint!(set, maximum(temp_cell_facet_neighbors_length)*(2^(Dim-1)))
    empty!(topology.neighbors)
    for (i, cell) in enumerate(grid.kopp_cells)
        new_i = refinement_cache.old_cell_to_new_cell_map[i]
        new_i == 0 && continue
        if refinement_cache.marked_for_coarsening[new_i]
            # TODO: Avoid using sets
            for facet in 1:NFacets
                for new_seq in 1:2^Dim
                    if ref_cell_neighborhood[new_seq][facet] == 0
                        child_neighbors_offset = temp_cell_facet_neighbors_offset[facet, i+new_seq]
                        child_neighbors_length = temp_cell_facet_neighbors_length[facet, i+new_seq]
                        child_neighbors = @view topology.neighbors_prev[child_neighbors_offset : child_neighbors_offset + child_neighbors_length - 1]
                        for child_neighbor in child_neighbors
                            parent_idx = grid.kopp_cells[child_neighbor[1]].parent
                            parent_idx_new = parent_idx < 0 ? parent_idx : refinement_cache.old_cell_to_new_cell_map[parent_idx]
                            if parent_idx < 0 || (refinement_cache.old_cell_to_new_cell_map[child_neighbor[1]] != 0 && !refinement_cache.marked_for_coarsening[parent_idx_new])
                                push!(set, FacetIndex(refinement_cache.old_cell_to_new_cell_map[child_neighbor[1]], child_neighbor[2]))
                            else
                                push!(set, FacetIndex(parent_idx_new, child_neighbor[2]))
                            end
                        end
                    end
                end
                topology.cell_facet_neighbors_offset[facet, new_i] = length(set) == 0 ? 0 : n_neighborhoods + 1
                topology.cell_facet_neighbors_length[facet, new_i] = length(set) == 0 ? 0 : length(set)
                n_neighborhoods += length(set)
                append!(topology.neighbors, set)
                empty!(set)
            end
        else
            for facet in 1:NFacets
                child_neighbors_offset = temp_cell_facet_neighbors_offset[facet, i]
                child_neighbors_length = temp_cell_facet_neighbors_length[facet, i]
                child_neighbors = @view topology.neighbors_prev[child_neighbors_offset : child_neighbors_offset + child_neighbors_length - 1]
                for child_neighbor in child_neighbors
                    parent_idx = grid.kopp_cells[child_neighbor[1]].parent
                    parent_idx_new = parent_idx < 0 ? parent_idx : refinement_cache.old_cell_to_new_cell_map[parent_idx]
                    if parent_idx < 0 || (refinement_cache.old_cell_to_new_cell_map[parent_idx] != 0 && refinement_cache.old_cell_to_new_cell_map[child_neighbor[1]] != 0 && !refinement_cache.marked_for_coarsening[parent_idx_new])
                        push!(set, FacetIndex(refinement_cache.old_cell_to_new_cell_map[child_neighbor[1]], child_neighbor[2]))
                    else
                        push!(set, FacetIndex(parent_idx_new, child_neighbor[2]))
                    end
                end
                topology.cell_facet_neighbors_offset[facet, new_i] = length(set) == 0 ? 0 : n_neighborhoods + 1
                topology.cell_facet_neighbors_length[facet, new_i] = length(set) == 0 ? 0 : length(set)
                n_neighborhoods += length(set)
                append!(topology.neighbors, set)
                empty!(set)
            end
        end
    end
    return n_neighborhoods
end

function update_coarsened_cells!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache) where Dim
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old cells
        new == 0 && continue
        old_cell = grid.kopp_cells_prev[old]
        grid.kopp_cells[new] =  KoppCell{Dim, Int}((old_cell.parent > 0 ? refinement_cache.old_cell_to_new_cell_map[old_cell.parent] : old_cell.parent) , old_cell.sequence, refinement_cache.marked_for_coarsening[new] ? true : old_cell.isleaf)
    end
end

function update_coarsened_koppcache!(
    grid::KoppGrid,
    refinement_cache::KoppRefinementCache,
    topology::KoppTopology,
    kopp_cache::AbstractAMRSynchronizer,
    kopp_values::ValuesCache,
    dh::DofHandler,
    temp_dh::DofHandler,
    NFacets::Int,
    Dim::Int)
    @inbounds for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old cells
        new == 0 && continue
        kopp_cache.cell_matrices[:,:,new] .= @view kopp_cache.cell_matrices[:,:,old]
        old != new && (kopp_cache.cell_matrices[:,:,new] .= 0.0)
    end
    @inbounds for cell_cache in Ferrite.CellIterator(grid)
        cell_idx = cell_cache.cellid
        cell = getcells(grid, cell_idx)
        if all(x -> x == 0.0, @view kopp_cache.cell_matrices[:,:,cell_idx])
            reinit!(kopp_values.cell_values, cell_cache)
            assemble_element_matrix!((@view kopp_cache.cell_matrices[:,:,cell_idx]), kopp_values)
        end
    end
    ninterfaces = 1
    facet_idx = 1
    nneighbors = 1
    interface_index = 1
    ii = Ferrite.InterfaceIterator(dh, topology)
    for (facet_idx, nneighbors, ninterfaces) in ii
        interface_cache = ii.cache
        _facet_idx = nneighbors == 1 ? facet_idx - 1 : facet_idx
        cell_idx = (_facet_idx - 1) ÷ (2*Dim) + 1
        facet_a = (_facet_idx - 1) % (2*Dim) + 1
        neighbor = FacetIndex(interface_cache.b.cc.cellid, interface_cache.b.current_facet_id)
        kopp_cache.interface_matrix_index[topology.cell_facet_neighbors_offset[facet_a, cell_idx] + nneighbors - 1] = ninterfaces
        rng  = topology.cell_facet_neighbors_offset[neighbor[2], neighbor[1]] : topology.cell_facet_neighbors_offset[neighbor[2], neighbor[1]] + topology.cell_facet_neighbors_length[neighbor[2], neighbor[1]] - 1
        # @assert any(x -> x == 0, @view kopp_cache.interface_matrix_index[rng])
        for j in rng
            kopp_cache.interface_matrix_index[j] != 0 && continue
            kopp_cache.interface_matrix_index[j] = ninterfaces
            break
        end
        reinit!(kopp_values.interface_values, interface_cache, topology)
        assemble_interface_matrix!((@view kopp_cache.interface_matrices[:,:,interface_index]), kopp_values)
        interface_index += 1
    end
end
