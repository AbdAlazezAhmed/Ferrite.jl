function _resize_marked_for_refinement!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache,
    cellset::AbstractSet{CellIndex}) where {Dim}
    n_refined_cells = length(cellset)
    new_length = length(grid.kopp_cells) + (2^Dim) * n_refined_cells
    resize!(refinement_cache.marked_for_refinement, new_length)
    resize!(refinement_cache.marked_for_coarsening, new_length)
    resize!(refinement_cache.new_cell_to_old_cell_map, new_length)

    (@view refinement_cache.marked_for_refinement[length(grid.kopp_cells) + 1 : end]) .= false
    (@view refinement_cache.marked_for_coarsening[length(grid.kopp_cells) + 1 : end]) .= false
    refinement_cache.new_cell_to_old_cell_map .= 0
    return nothing
end

function _update_refinement_cache_isactive!(
    grid::KoppGrid{Dim},
    topology::KoppTopology, # Old topology
    refinement_cache::KoppRefinementCache,
    cellset::AbstractSet{CellIndex}
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
                    get_refinement_level(grid.kopp_cells[neighbor.idx[1]]) > get_refinement_level(cell) && continue
                    # TODO: Invalidate interface matrices
                    # grid.interfaces_recompute[interface:interface + 1] .= true
                    n_refined_interfaces += 1
                end
            end
        end
        refinement_cache.new_cell_to_old_cell_map[refinement_cache.old_cell_to_new_cell_map[i]] = i
    end
    # for (i, cell) in enumerate(grid.kopp_cells)
    #     refinement_cache.new_cell_to_old_cell_map[refinement_cache.old_cell_to_new_cell_map[i]] = i
    # end
    old_length = length(refinement_cache.interfaces_data_updated_indices)
    # resize!(refinement_cache.interfaces_data_updated_indices, old_length + n_refined_interfaces * ((2^(Dim-1) - 1)) + 2Dim * length(cellset))
    # refinement_cache.interfaces_data_updated_indices .= 0
    old_length = length(refinement_cache.interfaces_updated_indices)

    return nothing
end
function add_inherited_neighbors!(n_neighborhoods::Int, cell_idx::Int, refinement_cache::KoppRefinementCache, facet::Int, topology::KoppTopology, grid::KoppGrid{Dim}, new_seq::Int) where {Dim}
    neighbors_count = 0
    parent_facet_neighborhood_offset = topology.cell_facet_neighbors_offset_prev[facet, cell_idx]
    parent_facet_neighborhood_length = topology.cell_facet_neighbors_length_prev[facet, cell_idx]
    parent_facet_neighbors = @view topology.neighbors_prev[parent_facet_neighborhood_offset:parent_facet_neighborhood_offset+parent_facet_neighborhood_length-1]
    for neighbor in parent_facet_neighbors
        neighbor_idx = neighbor[1]
        neighbor_facet = neighbor[2]
        cell = grid.kopp_cells_prev[cell_idx]
        neighbor_cell = grid.kopp_cells_prev[neighbor_idx]
        cell_facet_orientation_info = topology.root_facet_orientation_info[topology.root_idx_prev[cell_idx]+facet-1]
        neighbor_facet_orientation_info = topology.root_facet_orientation_info[topology.root_idx_prev[neighbor_idx]+neighbor_facet-1]
        cell_children = Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[facet]
        neighbor_children = SVector(Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[neighbor_facet])
        flip_interface = (cell_facet_orientation_info.flipped ⊻ neighbor_facet_orientation_info.flipped)
        shift = cell_facet_orientation_info.shift_index - (flip_interface ? neighbor_facet_orientation_info.shift_index : length(neighbor_children) - neighbor_facet_orientation_info.shift_index)
        new_neighbor_children = neighbor_children
        # URGERNT TODO: find a way to enable next line without allocations
        new_neighbor_children = reverse(neighbor_children)
        # new_neighbor_children = circshift((flip_interface ? reverse(neighbor_children) : neighbor_children), shift)
        # neighbor_local_seq = neighbor_cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * get_refinement_level(cell)))
        # neighbor_local_seq = (neighbor_cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * (max(0, get_refinement_level(cell)) )))) >>  ((Dim + 1) * (max(0, get_refinement_level(cell))))

        if (get_refinement_level(neighbor_cell) > get_refinement_level(cell))
            neighbor_local_seq = (neighbor_cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * (get_refinement_level(neighbor_cell) - get_refinement_level(cell) -1)))) >>  ((Dim + 1) * (get_refinement_level(neighbor_cell) - get_refinement_level(cell) -1))
            if findfirst(==(new_seq), cell_children) == findfirst(==(neighbor_local_seq), new_neighbor_children)
                if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor_idx]]
                    topology.neighbors[n_neighborhoods + neighbors_count + 1 : n_neighborhoods + neighbors_count + length(neighbor_children)] = [FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor_idx] + neighbor_child, neighbor_facet) for neighbor_child in neighbor_children]
                    neighbors_count += length(neighbor_children)
                    continue
                else
                    topology.neighbors[n_neighborhoods + neighbors_count + 1] = FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor_idx], neighbor_facet)
                    neighbors_count += 1
                    continue
                end
            end
        elseif (get_refinement_level(neighbor_cell) < get_refinement_level(cell))
            _new_seq = (cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * (get_refinement_level(cell) - get_refinement_level(neighbor_cell) -1)))) >>  ((Dim + 1) * ((get_refinement_level(cell) - get_refinement_level(neighbor_cell) -1)))
            # if findfirst(==(_new_seq), cell_children) == findfirst(==(neighbor_local_seq), new_neighbor_children)
                if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor_idx]]
                    topology.neighbors[n_neighborhoods + neighbors_count + 1] = FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor_idx] + new_neighbor_children[findfirst(==(_new_seq), cell_children)], neighbor_facet)
                    neighbors_count += 1
                    continue
                else
                    topology.neighbors[n_neighborhoods + neighbors_count + 1] = FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor_idx], neighbor_facet)
                    neighbors_count += 1
                    continue
                end
            # end
        elseif (get_refinement_level(neighbor_cell) == get_refinement_level(cell))
            # _new_seq = (cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * (max(get_refinement_level(neighbor_cell) - 1, 0))))) >>  ((Dim + 1) * (max(get_refinement_level(neighbor_cell), 1) - 1))
            # neighbor_cell.sequence == 0 &&
            if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor_idx]]
                topology.neighbors[n_neighborhoods + neighbors_count + 1] = FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor_idx] + new_neighbor_children[findfirst(==(new_seq), cell_children)], neighbor_facet)
                neighbors_count += 1
                continue
            else
                topology.neighbors[n_neighborhoods + neighbors_count + 1] = FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor_idx], neighbor_facet)
                neighbors_count += 1
                continue
            end
        end
    end

    return neighbors_count
end
function update_neighbors!(
    grid::KoppGrid{Dim},
    topology::KoppTopology,
    refinement_cache::KoppRefinementCache) where {Dim}
    NFacets = 2 * Dim
    temp_cell_facet_neighbors_offset = topology.cell_facet_neighbors_offset_prev
    temp_cell_facet_neighbors_length = topology.cell_facet_neighbors_length_prev
    n_neighborhoods = 0
    ref_cell_neighborhood = get_ref_cell_neighborhood(getrefshape(getcelltype(grid.base_grid)))
    for (i, cell) in enumerate(grid.kopp_cells_prev)
        new_i = refinement_cache.old_cell_to_new_cell_map[i]
        if refinement_cache.marked_for_refinement[new_i]
            for new_seq in 1:2^Dim
                for facet in 1:NFacets
                    n_facet_neighbors = 0
                    if ref_cell_neighborhood[new_seq][facet] == 0
                        n_inherited_neighbors = add_inherited_neighbors!(n_neighborhoods, i, refinement_cache, facet, topology, grid, new_seq)
                        n_facet_neighbors += n_inherited_neighbors
                    else
                        n_facet_neighbors += 1
                        topology.neighbors[n_neighborhoods + 1] = FacetIndex(new_i + ref_cell_neighborhood[new_seq][facet], get_ref_face_neighborhood(getrefshape(getcelltype(grid.base_grid)))[facet])
                    end
                    n_neighborhoods += n_facet_neighbors
                end
            end
        else
            for facet in 1:NFacets
                n_facet_neighbors = 0
                neighbors_offset = temp_cell_facet_neighbors_offset[facet, i]
                neighbors_length = temp_cell_facet_neighbors_length[facet, i]
                cell_children = Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[facet]

                # BUG SOLUTION HERE: make it into cases
                for neighbor in @view topology.neighbors_prev[neighbors_offset:neighbors_offset+neighbors_length-1]
                    neighbor_facet = neighbor[2]
                    neighbor_idx = neighbor[1]
                    neighbor_children = SVector(Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[neighbor_facet])
                    cell_facet_orientation_info = topology.root_facet_orientation_info[topology.root_idx_prev[i]+facet-1]
                    neighbor_facet_orientation_info = topology.root_facet_orientation_info[topology.root_idx_prev[neighbor[1]]+neighbor_facet-1]
                    flip_interface = (cell_facet_orientation_info.flipped ⊻ neighbor_facet_orientation_info.flipped)
                    shift = cell_facet_orientation_info.shift_index - (flip_interface ? neighbor_facet_orientation_info.shift_index : length(neighbor_children) - neighbor_facet_orientation_info.shift_index)
                    new_neighbor_children = neighbor_children
                    # URGERNT TODO: find a way to enable next line without allocations
                    new_neighbor_children = reverse(neighbor_children)
                    if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor[1]]]
                        neighbor_cell = grid.kopp_cells_prev[neighbor[1]]
                        if get_refinement_level(neighbor_cell) < get_refinement_level(cell)
                            _new_seq = (cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * (get_refinement_level(cell) - get_refinement_level(neighbor_cell) -1)))) >>  ((Dim + 1) * ((get_refinement_level(cell) - get_refinement_level(neighbor_cell) -1)))
                            topology.neighbors[n_neighborhoods + n_facet_neighbors + 1] = FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor_idx] + new_neighbor_children[findfirst(==(_new_seq), cell_children)], neighbor_facet)
                            n_facet_neighbors += 1
                        else
                            # TODO: double check if we use old or new index here
                            for (i, new_seq) in enumerate(Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[neighbor[2]])
                                topology.neighbors[n_neighborhoods + n_facet_neighbors + i] = FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor[1]] + new_seq, neighbor[2])
                            end
                            n_facet_neighbors += 2^(Dim - 1)
                        end


                    else
                        n_facet_neighbors += 1
                        topology.neighbors[n_neighborhoods + n_facet_neighbors] = FacetIndex(refinement_cache.old_cell_to_new_cell_map[neighbor[1]], neighbor[2])
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

function update_cells!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache) where Dim
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old cells
        old_cell = grid.kopp_cells_prev[old]
        grid.kopp_cells[new] =  KoppCell{Dim, Int}((old_cell.parent > 0 ? refinement_cache.old_cell_to_new_cell_map[old_cell.parent] : old_cell.parent) , old_cell.sequence, old_cell.isleaf)
    end
    for (cell_idx, cell) in enumerate(grid.kopp_cells)
        refinement_cache.marked_for_refinement[cell_idx] || continue
        grid.kopp_cells[cell_idx] = KoppCell{Dim, Int}(cell.parent, cell.sequence, false)
        for new_seq in 1 : 2^Dim
            child_idx = cell_idx + new_seq
            grid.kopp_cells[child_idx] = KoppCell{Dim, Int}(cell_idx, (cell.sequence << (Dim + 1)) + new_seq, true)
        end
    end
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
    for (i, cell) in enumerate(grid.kopp_cells)
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
                cell_children = Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[facet]
                for neighbor in @view topology.neighbors_prev[neighbors_offset:neighbors_offset+neighbors_length-1]
                    neighbor_facet = neighbor[2]
                    neighbor_idx = neighbor[1]
                    neighbor_children = SVector(Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[neighbor_facet])
                    cell_facet_orientation_info = topology.root_facet_orientation_info[topology.root_idx_prev[i]+facet-1]
                    neighbor_facet_orientation_info = topology.root_facet_orientation_info[topology.root_idx_prev[neighbor[1]]+neighbor_facet-1]
                    flip_interface = (cell_facet_orientation_info.flipped ⊻ neighbor_facet_orientation_info.flipped)
                    shift = cell_facet_orientation_info.shift_index - (flip_interface ? neighbor_facet_orientation_info.shift_index : length(neighbor_children) - neighbor_facet_orientation_info.shift_index)
                    new_neighbor_children = neighbor_children
                    # URGERNT TODO: find a way to enable next line without allocations
                    new_neighbor_children = reverse(neighbor_children)
                    if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor[1]]]
                        neighbor_cell = grid.kopp_cells_prev[neighbor[1]]
                        if get_refinement_level(neighbor_cell) < get_refinement_level(cell)
                            n_facet_neighbors += 1
                        else
                            n_facet_neighbors += 2^(Dim - 1)
                        end
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
reference_coordinates_static(::Val{2}, ::Val{2}) = (nothing,)
reference_coordinates_static(::Val{2}, ::Val{1}) = (Vec{2, Float64}((-1.0, -1.0)),
Vec{2, Float64}((1.0, -1.0)),
Vec{2, Float64}((1.0, 1.0)),
Vec{2, Float64}((-1.0, 1.0)),)
function get_ref_inherited_dofs!(children_dofs, parent_dofs, ::Val{Dim}, ::Val{Order}) where {Dim, Order}
    coords = reference_coordinates_static(Val(Dim), Val(Order))
    root_coords = ntuple(i -> Node(coords[i]), length(coords))
    dofs = zero(MVector{(Order+1)^Dim * 2^(Dim), Int})
    j = 0
    @inbounds for i in 1:2^Dim
        child_coords = _process_seq(Dim, root_coords, i)
        for dof in child_coords
            j+=1
            x = findfirst(x -> x == dof, root_coords)
            isnothing(x) && continue
            dofs[j] = x
        end
    end
    new_dofs = maximum(parent_dofs)
    @inbounds for i in eachindex(children_dofs)
        if iszero(dofs[i])
            children_dofs[i] =  new_dofs + 1
            new_dofs += 1
        else
            children_dofs[i] = parent_dofs[dofs[i]]
        end
    end
    return nothing
end

function update_dofs_H1!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache,
    dh::DofHandler,
    temp_celldofs::Vector{Int},
    temp_cell_dofs_offset::Vector{Int},
    temp_cell_to_subdofhandler::Vector{Int},
    ::Val{Order}) where {Dim, Order}
    dh.cell_dofs_offset .= 0
    # dh.cell_dofs .= 1
    @inbounds for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        new == 0 && continue
        dh.cell_dofs_offset[new] = temp_cell_dofs_offset[old]
    end
    refined = 0
    coarsened = 0
    @inbounds @views for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        new == 0 && continue
        ndofs_per_cell = dh.subdofhandlers[temp_cell_to_subdofhandler[old]].ndofs_per_cell
        dh.cell_dofs_offset[new] += (2^Dim) * ndofs_per_cell * refined
        dh.cell_dofs_offset[new] -= (2^Dim) * ndofs_per_cell * coarsened
        if refinement_cache.marked_for_refinement[new]
            # for i in new + 1 : length(dh.cell_dofs_offset)
                # dh.cell_dofs_offset[i] += (2^Dim) * ndofs_per_cell
            # end
            refined += 1
        end
        if refinement_cache.marked_for_coarsening[new]
            coarsened += 1
        end
        dh.cell_dofs[dh.cell_dofs_offset[new] : dh.cell_dofs_offset[new] + dh.subdofhandlers[1].ndofs_per_cell - 1] .= temp_celldofs[temp_cell_dofs_offset[old] : temp_cell_dofs_offset[old] + dh.subdofhandlers[1].ndofs_per_cell - 1]
        dh.cell_to_subdofhandler[new] = temp_cell_to_subdofhandler[old]
    end
    max_dof = 0
    refined = 0
    @inbounds @views for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        new == 0 && continue
        # copy old dofhander vectors
        ndofs_per_cell = dh.subdofhandlers[temp_cell_to_subdofhandler[old]].ndofs_per_cell
        if refinement_cache.marked_for_refinement[new]
            cell_dofs_current_cell = dh.cell_dofs[dh.cell_dofs_offset[new] :  dh.cell_dofs_offset[new] + dh.subdofhandlers[1].ndofs_per_cell - 1]
            max_dof = maximum(cell_dofs_current_cell)
            for i in dh.cell_dofs_offset[new] + (2^Dim + 1) * ndofs_per_cell : length(dh.cell_dofs)
                dh.cell_dofs[i] > max_dof || continue
                dh.cell_dofs[i] += (2^Dim - 1) * ndofs_per_cell
            end
        end
        if refinement_cache.marked_for_coarsening[new]
            cell_dofs_current_cell = dh.cell_dofs[dh.cell_dofs_offset[new] :  dh.cell_dofs_offset[new] + dh.subdofhandlers[1].ndofs_per_cell - 1]
            max_dof = maximum(cell_dofs_current_cell)
            for i in dh.cell_dofs_offset[new] : length(dh.cell_dofs)
                dh.cell_dofs[i] > max_dof || continue
                dh.cell_dofs[i] -= (2^Dim - 1) * ndofs_per_cell
            end
        end
    end

    @views for (old, cell_idx) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        cell_idx == 0 && continue
        # ASSUMPTION: single field
        if refinement_cache.marked_for_refinement[cell_idx]
            ndofs_per_cell = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx]].ndofs_per_cell
            @inbounds for new_seq in 1 : 2^Dim
                child_idx = cell_idx + new_seq
                # DofHandler add new dofs and offsets
                # For DG only
                dh.cell_dofs_offset[child_idx] = dh.cell_dofs_offset[cell_idx] + new_seq * ndofs_per_cell
                dh.cell_to_subdofhandler[child_idx] = dh.cell_to_subdofhandler[cell_idx]
                dh.ndofs += dh.subdofhandlers[1].ndofs_per_cell
            end
            children_dofs = dh.cell_dofs[dh.cell_dofs_offset[cell_idx+1] : dh.cell_dofs_offset[cell_idx+1] + (2^Dim) * ndofs_per_cell - 1]
            parent_dofs = dh.cell_dofs[dh.cell_dofs_offset[cell_idx] : dh.cell_dofs_offset[cell_idx + 1] - 1]

            get_ref_inherited_dofs!(children_dofs, parent_dofs, Val(Dim), Val(Order))
        elseif refinement_cache.marked_for_coarsening[cell_idx]
        end

    end

end


function update_dofs_L2!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache,
    dh::DofHandler,
    temp_celldofs::Vector{Int},
    temp_cell_dofs_offset::Vector{Int},
    temp_cell_to_subdofhandler::Vector{Int},
    ::Val{Order}) where {Dim, Order}
    dh.cell_dofs_offset .= 0
    # dh.cell_dofs .= 1
    @inbounds for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        new == 0 && continue
        dh.cell_dofs_offset[new] = temp_cell_dofs_offset[old]
    end
    refined = 0
    coarsened = 0
    @inbounds @views for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        new == 0 && continue
        ndofs_per_cell = dh.subdofhandlers[temp_cell_to_subdofhandler[old]].ndofs_per_cell
        dh.cell_dofs_offset[new] += (2^Dim) * ndofs_per_cell * refined
        dh.cell_dofs_offset[new] -= (2^Dim) * ndofs_per_cell * coarsened
        if refinement_cache.marked_for_refinement[new]
            # for i in new + 1 : length(dh.cell_dofs_offset)
                # dh.cell_dofs_offset[i] += (2^Dim) * ndofs_per_cell
            # end
            refined += 1
        end
        if refinement_cache.marked_for_coarsening[new]
            coarsened += 1
        end
        dh.cell_dofs[dh.cell_dofs_offset[new] : dh.cell_dofs_offset[new] + dh.subdofhandlers[1].ndofs_per_cell - 1] .= temp_celldofs[temp_cell_dofs_offset[old] : temp_cell_dofs_offset[old] + dh.subdofhandlers[1].ndofs_per_cell - 1]
        dh.cell_to_subdofhandler[new] = temp_cell_to_subdofhandler[old]
    end

    @views for (old, cell_idx) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        cell_idx == 0 && continue
        # ASSUMPTION: single field
        if refinement_cache.marked_for_refinement[cell_idx]
            ndofs_per_cell = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx]].ndofs_per_cell
            @inbounds for new_seq in 1 : 2^Dim
                child_idx = cell_idx + new_seq
                # DofHandler add new dofs and offsets
                # For DG only
                dh.cell_dofs_offset[child_idx] = dh.cell_dofs_offset[cell_idx] + new_seq * ndofs_per_cell
                dh.cell_to_subdofhandler[child_idx] = dh.cell_to_subdofhandler[cell_idx]
                dh.ndofs += dh.subdofhandlers[1].ndofs_per_cell
            end
            # children_dofs = dh.cell_dofs[dh.cell_dofs_offset[cell_idx+1] : dh.cell_dofs_offset[cell_idx+1] + (2^Dim) * ndofs_per_cell - 1]
            # parent_dofs = dh.cell_dofs[dh.cell_dofs_offset[cell_idx] : dh.cell_dofs_offset[cell_idx + 1] - 1]

            # get_ref_inherited_dofs!(children_dofs, parent_dofs, Val(Dim), Val(Order))
        elseif refinement_cache.marked_for_coarsening[cell_idx]
        end

    end

    max_dof = 0
    current_dof = 1
    dh.cell_dofs .= 0
    # TODO: URGENT HARDCODED
    ip = Lagrange{RefQuadrilateral, 1}()
    @views for cc in CellIterator(grid)
        cell_idx = cellid(cc)
        cell = grid.kopp_cells[cell_idx]
        cell.isleaf || continue
        facet_dofs_idx = Ferrite.facetdof_indices(ip)
        cell_dofs = dh.cell_dofs[dh.cell_dofs_offset[cell_idx] : dh.cell_dofs_offset[cell_idx] + dh.subdofhandlers[1].ndofs_per_cell - 1]

        for facet in 1:2*Dim
            neighborhood = getneighborhood(topology, FacetIndex(cell_idx, facet))
            length(neighborhood) > 1 && continue
            for (i, dof_idx) in enumerate(facet_dofs_idx[facet])
                if iszero(cell_dofs[dof_idx])
                    cell_dofs[dof_idx] = current_dof
                    if length(neighborhood) == 1
                        neighbor = neighborhood[]
                        # neighbor[1] < cell_idx && continue
                        get_refinement_level(grid.kopp_cells[neighbor[1]]) < get_refinement_level(cell) && continue
                        neighbor_dofs = dh.cell_dofs[dh.cell_dofs_offset[neighbor[1]] : dh.cell_dofs_offset[neighbor[1]] + dh.subdofhandlers[1].ndofs_per_cell - 1]
                        #TODO URGENT 3D take care of the transformation
                        neighbor_facet_dofs = facet_dofs_idx[neighbor[2]]
                        iszero(neighbor_dofs[neighbor_facet_dofs[end - i + 1]]) || continue
                        neighbor_dofs[neighbor_facet_dofs[end - i + 1]] = current_dof
                    end
                    current_dof += 1
                else
                    if length(neighborhood) == 1
                        neighbor = neighborhood[]
                        neighbor[1] < cell_idx && continue
                        get_refinement_level(grid.kopp_cells[neighbor[1]]) < get_refinement_level(cell) && continue
                        neighbor_dofs = dh.cell_dofs[dh.cell_dofs_offset[neighbor[1]] : dh.cell_dofs_offset[neighbor[1]] + dh.subdofhandlers[1].ndofs_per_cell - 1]
                        #TODO URGENT 3D take care of the transformation
                        neighbor_facet_dofs = facet_dofs_idx[neighbor[2]]
                        neighbor_dofs[neighbor_facet_dofs[end - i + 1]] = current_dof
                    end
                end
            end
        end
        parent_cell_idx = cell.parent
        for i in 1:get_refinement_level(cell)
            parent = grid.kopp_cells[parent_cell_idx]
            cell_seq = (cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * ( i -1)))) >>  ((Dim + 1) * (i -1))
            parent_dofs = dh.cell_dofs[dh.cell_dofs_offset[parent_cell_idx] : dh.cell_dofs_offset[parent_cell_idx] + dh.subdofhandlers[1].ndofs_per_cell - 1]
            iszero(parent_dofs[cell_seq]) && (parent_dofs[cell_seq] = cell_dofs[cell_seq])
            for facet in 1:2*Dim
                neighborhood = getneighborhood(topology, FacetIndex(parent_cell_idx, facet))
                length(neighborhood) > 1 && continue
                for (i, dof_idx) in enumerate(facet_dofs_idx[facet])
                    dof_idx != cell_seq && continue
                    if length(neighborhood) == 1
                        neighbor = neighborhood[]
                        get_refinement_level(grid.kopp_cells[neighbor[1]]) < get_refinement_level(parent) && continue
                        neighbor_dofs = dh.cell_dofs[dh.cell_dofs_offset[neighbor[1]] : dh.cell_dofs_offset[neighbor[1]] + dh.subdofhandlers[1].ndofs_per_cell - 1]
                        #TODO URGENT 3D take care of the transformation
                        neighbor_facet_dofs = facet_dofs_idx[neighbor[2]]
                        iszero(neighbor_dofs[neighbor_facet_dofs[end - i + 1]]) && (neighbor_dofs[neighbor_facet_dofs[end - i + 1]] = current_dof)
                    end
                end
            end
            parent_seq = (cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * ( i+1 -1)))) >>  ((Dim + 1) * (i+1 -1))
            (parent_seq != 0 || parent_seq != cell_seq) && break
            parent_cell_idx = parent.parent
        end
    end
    # @inbounds @views for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
    #     new == 0 && continue
    #     # copy old dofhander vectors
    #     ndofs_per_cell = dh.subdofhandlers[temp_cell_to_subdofhandler[old]].ndofs_per_cell
    #     if refinement_cache.marked_for_refinement[new]
    #         cell_dofs_current_cell = dh.cell_dofs[dh.cell_dofs_offset[new] :  dh.cell_dofs_offset[new] + dh.subdofhandlers[1].ndofs_per_cell - 1]
    #         max_dof = maximum(cell_dofs_current_cell)
    #         for i in dh.cell_dofs_offset[new] + (2^Dim + 1) * ndofs_per_cell : length(dh.cell_dofs)
    #             dh.cell_dofs[i] > max_dof || continue
    #             dh.cell_dofs[i] += (2^Dim - 1) * ndofs_per_cell
    #         end
    #     end
    #     if refinement_cache.marked_for_coarsening[new]
    #         cell_dofs_current_cell = dh.cell_dofs[dh.cell_dofs_offset[new] :  dh.cell_dofs_offset[new] + dh.subdofhandlers[1].ndofs_per_cell - 1]
    #         max_dof = maximum(cell_dofs_current_cell)
    #         for i in dh.cell_dofs_offset[new] : length(dh.cell_dofs)
    #             dh.cell_dofs[i] > max_dof || continue
    #             dh.cell_dofs[i] -= (2^Dim - 1) * ndofs_per_cell
    #         end
    #     end
    # end



end
