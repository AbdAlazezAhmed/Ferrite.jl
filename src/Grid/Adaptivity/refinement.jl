function _resize_marked_for_refinement!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache,
    cellset::Set{CellIndex}) where {Dim}
    n_refined_cells = length(cellset)
    new_length = length(grid.kopp_cells) + (2^Dim) * n_refined_cells
    resize!(refinement_cache.marked_for_refinement, new_length)
    resize!(refinement_cache.new_cell_to_old_cell_map, new_length)

    refinement_cache.marked_for_refinement[length(grid.kopp_cells) + 1 : end] .= false
    refinement_cache.new_cell_to_old_cell_map .= 0
    return nothing
end

function _update_refinement_cache_isactive!(
    grid::KoppGrid{Dim},
    topology::KoppTopology, # Old topology
    refinement_cache::KoppRefinementCache,
    cellset::Set{CellIndex}
    ) where {Dim}
    n_refined_interfaces = 0
    for (i, cell) in enumerate(grid.kopp_cells)
        if CellIndex(i) ∈ cellset
            refinement_cache.old_cell_to_new_cell_map[i+1:end] .= (@view refinement_cache.old_cell_to_new_cell_map[i+1:end]) .+ 2^Dim
            refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[i]] = true
            # kopp_cache.ansatz_isactive[@view dh.cell_dofs[dh.cell_dofs_offset[i]:dh.cell_dofs_offset[i]+dh.subdofhandlers[1].ndofs_per_cell-1]] .= false
            # kopp_cache.ansatz_isactive[dh.cell_dofs_offset[i]+dh.subdofhandlers[1].ndofs_per_cell : dh.cell_dofs_offset[i]+dh.subdofhandlers[1].ndofs_per_cell*(2^Dim) - 1] .= true
            for facet in 1:2Dim
                interface_offset = topology.cell_facet_neighbors_offset[facet, i]
                interface_length = topology.cell_facet_neighbors_length[facet, i]
                for interface_index_idx in interface_offset : interface_offset + interface_length - 1
                    interface_offset == 0 && continue
                    neighbor = topology.neighbors[interface_index_idx]
                    # Refine interface when the current cell is lower index than the neighbor only to avoid race conditions
                    get_refinement_level(grid.kopp_cells[neighbor.idx[1]]) == get_refinement_level(cell) && i > neighbor.idx[1] && continue
                    # Refine the interface only if the neighbor is not finer than the current cell
                    get_refinement_level(grid.kopp_cells[neighbor.idx[1]]) < get_refinement_level(cell) && continue
                    # TODO: Invalidate interface matrices
                    # grid.interfaces_recompute[interface:interface + 1] .= true
                    n_refined_interfaces += 1
                end
            end
        end
        refinement_cache.new_cell_to_old_cell_map[refinement_cache.old_cell_to_new_cell_map[i]] = i
    end
    old_length = length(refinement_cache.interfaces_data_updated_indices)
    resize!(refinement_cache.interfaces_data_updated_indices, old_length + n_refined_interfaces * ((2^(Dim-1) - 1)) + 2Dim * length(cellset))
    refinement_cache.interfaces_data_updated_indices .= 0
    old_length = length(refinement_cache.interfaces_updated_indices)
    resize!(refinement_cache.interfaces_updated_indices, old_length + n_refined_interfaces * ((2^(Dim-1) - 1)) + 2Dim * length(cellset))
    refinement_cache.interfaces_data_updated_indices .= 0
    return nothing
end

function count_neighbors_update_indexing!(
    grid::KoppGrid{Dim},
    topology::KoppTopology,
    refinement_cache::KoppRefinementCache) where {Dim}
    NFacets = 2 * Dim
    temp_cell_facet_neighbors_offset = topology.cell_facet_neighbors_offset_prev
    temp_cell_facet_neighbors_length = topology.cell_facet_neighbors_length_prev
    n_neighborhoods = 0
    ref_cell_neighborhood = get_ref_cell_neighborhood(getrefshape(getcelltype(grid.base_grid)))
    for i in 1:length(grid.kopp_cells)
        new_i = refinement_cache.old_cell_to_new_cell_map[i]
        if refinement_cache.marked_for_refinement[new_i]
            for new_seq in 1:2^Dim
                for facet in 1:NFacets
                    n_facet_neighbors = 0
                    if ref_cell_neighborhood[new_seq][facet] == 0
                        n_inherited_neighbors = count_inherited_neighbors(i, refinement_cache, facet, topology, grid, new_seq)
                        n_facet_neighbors += n_inherited_neighbors
                    else
                        n_facet_neighbors += 1
                    end
                    topology.cell_facet_neighbors_offset[facet, new_i+new_seq] = n_facet_neighbors == 0 ? 0 : n_neighborhoods + 1
                    topology.cell_facet_neighbors_length[facet, new_i+new_seq] = n_facet_neighbors == 0 ? 0 : n_facet_neighbors
                    n_neighborhoods += n_facet_neighbors
                end
            end
        else
            for facet in 1:NFacets
                n_facet_neighbors = 0
                neighbors_offset = temp_cell_facet_neighbors_offset[facet, i]
                neighbors_length = temp_cell_facet_neighbors_length[facet, i]
                for neighbor in @view topology.neighbors_prev[neighbors_offset:neighbors_offset+neighbors_length-1]
                    if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor[1]]]
                        n_facet_neighbors += 2^(Dim - 1)
                    else
                        n_facet_neighbors += 1
                    end
                end
                topology.cell_facet_neighbors_offset[facet, new_i] = n_facet_neighbors == 0 ? 0 : n_neighborhoods + 1
                topology.cell_facet_neighbors_length[facet, new_i] = n_facet_neighbors == 0 ? 0 : n_facet_neighbors
                n_neighborhoods += n_facet_neighbors
            end
        end
    end
    return n_neighborhoods
end
# TODO: URGENT OPTIMIZE
function get_ref_inherited_dofs!(children_dofs, parent_dofs, interpolation::Interpolation)
    coords = Ferrite.reference_coordinates(interpolation)
    Dim = Ferrite.getrefdim(interpolation)
    root_coords = ntuple(i -> Node(coords[i]), length(coords))
    dofs = zeros(Int, Ferrite.getnbasefunctions(interpolation) * 2^(Dim))
    j = 0
    for i in 1:2^Dim
        child_coords = _process_seq(Dim, root_coords, i)
        for dof in child_coords
            j+=1
            x = findfirst(x -> x == dof, root_coords)
            isnothing(x) && continue
            dofs[j] = x
        end
    end
    new_dofs = maximum(parent_dofs)
    for i in eachindex(children_dofs)
        if iszero(dofs[i])
            children_dofs[i] =  new_dofs + 1
            new_dofs += 1
        else
            children_dofs[i] = parent_dofs[dofs[i]]
        end
    end
    children_dofs
end

function update_dofs!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache,
    dh::DofHandler,
    temp_celldofs::Vector{Int},
    temp_cell_dofs_offset::Vector{Int},
    temp_cell_to_subdofhandler::Vector{Int}) where Dim
    dh.cell_dofs_offset .= 0
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        dh.cell_dofs_offset[new] = temp_cell_dofs_offset[old]
    end
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        ndofs_per_cell = dh.subdofhandlers[temp_cell_to_subdofhandler[old]].ndofs_per_cell
        if refinement_cache.marked_for_refinement[new]
            dh.cell_dofs_offset[new + 1 : end] .+= (2^Dim) * ndofs_per_cell
        end
        dh.cell_dofs[dh.cell_dofs_offset[new] : dh.cell_dofs_offset[new] + dh.subdofhandlers[1].ndofs_per_cell - 1] .= @view temp_celldofs[temp_cell_dofs_offset[old] : temp_cell_dofs_offset[old] + dh.subdofhandlers[1].ndofs_per_cell - 1]
        dh.cell_to_subdofhandler[new] = temp_cell_to_subdofhandler[old]
    end
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        ndofs_per_cell = dh.subdofhandlers[temp_cell_to_subdofhandler[old]].ndofs_per_cell
        if refinement_cache.marked_for_refinement[new]
            dh.cell_dofs[dh.cell_dofs_offset[new] + (2^Dim + 1) * ndofs_per_cell : end] .+= (2^Dim - 1) * ndofs_per_cell
        end
    end
    @info dh.cell_dofs
    @info dh.cell_dofs_offset
    @info refinement_cache.marked_for_refinement
    for (old, cell_idx) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        @info cell_idx
        refinement_cache.marked_for_refinement[cell_idx] || continue
        @info cell_idx

        ndofs_per_cell = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx]].ndofs_per_cell
        for new_seq in 1 : 2^Dim
            child_idx = cell_idx + new_seq
            # DofHandler add new dofs and offsets
            # For DG only
            dh.cell_dofs_offset[child_idx] = dh.cell_dofs_offset[cell_idx] + new_seq * ndofs_per_cell
            dh.cell_to_subdofhandler[child_idx] = dh.cell_to_subdofhandler[cell_idx]
            dh.ndofs += dh.subdofhandlers[1].ndofs_per_cell
        end
        children_dofs = @view dh.cell_dofs[dh.cell_dofs_offset[cell_idx+1] : dh.cell_dofs_offset[cell_idx+1] + 2^Dim * ndofs_per_cell - 1]
        parent_dofs = @view dh.cell_dofs[dh.cell_dofs_offset[cell_idx] : dh.cell_dofs_offset[cell_idx + 1] - 1]
        # ASSUMPTION: single field
        get_ref_inherited_dofs!(children_dofs, parent_dofs, dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx]].field_interpolations[1])
    end
    @info dh.cell_dofs
    @info dh.cell_dofs_offset
end
