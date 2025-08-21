
function _resize_topology!(topology::KoppTopology, new_length::Int, ::Val{Dim}) where {Dim}
    resize!(topology.cell_facet_neighbors_offset, 2Dim, new_length)
    resize!(topology.cell_facet_neighbors_length, 2Dim, new_length)
    return nothing
end

function zero_topology!(topology::KoppTopology)
    topology.cell_facet_neighbors_offset .= 0
    topology.cell_facet_neighbors_length .= 0
    return nothing
end

function _calc_interfaces_dict_prev(
    grid::KoppGrid{Dim},
    topology::KoppTopology) where {Dim}
    NFacets = 2 * Dim
    interfaces_dict_prev = Dict{FacetIndex, Int}()
    interfaces_dict_prev_vec = Vector{Pair{FacetIndex, Int}}(undef, length(topology.neighbors)÷2)
    # sizehint!(interfaces_dict_prev, length(topology.neighbors)÷2)
    i = 1
    # @time "iter" begin
        for ic in Ferrite.InterfaceIterator(grid, topology)
            interfaces_dict_prev_vec[i] = FacetIndex(Ferrite.cellid(ic.a), ic.a.current_facet_id) => i
            i += 1
        end
    # end
    interfaces_dict_prev = Dict(interfaces_dict_prev_vec)
    return interfaces_dict_prev
end

function _calc_interfaces_map(grid::KoppGrid{Dim},
    topology::KoppTopology,
    refinement_cache::KoppRefinementCache,
    dict_prev::Dict) where {Dim}
    interface_index = 1
    ii = Ferrite.InterfaceIterator(grid, topology)
    for ic in ii
        (refinement_cache.marked_for_refinement[Ferrite.cellid(ic.a)] ||
        refinement_cache.marked_for_refinement[Ferrite.cellid(ic.b)] ||
        refinement_cache.marked_for_coarsening[Ferrite.cellid(ic.a)] ||
        refinement_cache.marked_for_coarsening[Ferrite.cellid(ic.b)] ) && continue
        cell_idx = refinement_cache.new_cell_to_old_cell_map[Ferrite.cellid(ic.a)]
        iszero(cell_idx) && (interface_index += 1; continue)
        old_index = dict_prev[FacetIndex(cell_idx, ic.a.current_facet_id)]
        (old_index<=0) && (interface_index += 1; continue)
        refinement_cache.interfaces_data_updated_indices[old_index] = interface_index
        interface_index += 1
    end
end

function update_root_idx!(
    grid::KoppGrid{Dim},
    topology::KoppTopology,
    refinement_cache::KoppRefinementCache) where Dim
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old cells
        topology.root_idx[new] =  topology.root_idx_prev[old]
    end
    for (cell_idx, cell) in enumerate(grid.kopp_cells)
        refinement_cache.marked_for_refinement[cell_idx] || continue
        for new_seq in 1 : 2^Dim
            child_idx = cell_idx + new_seq
            topology.root_idx[child_idx] = topology.root_idx[cell_idx]
        end
    end
end

function update_coarsened_root_idx!(
    grid::KoppGrid{Dim},
    topology::KoppTopology,
    refinement_cache::KoppRefinementCache) where Dim
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old cells
        new == 0 && continue
        topology.root_idx[new] =  topology.root_idx_prev[old]
    end
end

function _resize_bunch_of_stuff!(
    grid::KoppGrid{Dim},
    topology::KoppTopology,
    n_neighborhoods::Int,
    new_length::Int) where {Dim}
    resize!(topology.neighbors, n_neighborhoods)
    resize!(grid.kopp_cells, new_length)
    resize!(topology.root_idx, new_length)
    return nothing
end

function count_inherited_neighbors(cell_idx::Int, refinement_cache::KoppRefinementCache, facet::Int, topology::KoppTopology, grid::KoppGrid{Dim}, new_seq::Int) where {Dim}
    neighbors_count = 0
    parent_facet_neighborhood_offset = topology.cell_facet_neighbors_offset_prev[facet, cell_idx]
    parent_facet_neighborhood_length = topology.cell_facet_neighbors_length_prev[facet, cell_idx]
    parent_facet_neighbors = @view topology.neighbors_prev[parent_facet_neighborhood_offset:parent_facet_neighborhood_offset+parent_facet_neighborhood_length-1]
    for neighbor in parent_facet_neighbors
        neighbor_idx = neighbor[1]
        neighbor_facet = neighbor[2]
        cell = grid.kopp_cells[cell_idx]
        neighbor_cell = grid.kopp_cells[neighbor_idx]
        cell_facet_orientation_info = topology.root_facet_orientation_info[topology.root_idx_prev[cell_idx]+facet-1]
        neighbor_facet_orientation_info = topology.root_facet_orientation_info[topology.root_idx_prev[neighbor_idx]+neighbor_facet-1]
        cell_children = Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[facet]
        neighbor_children = SVector(Ferrite.reference_facets(getrefshape(getcelltype(grid.base_grid)))[neighbor_facet])
        flip_interface = (cell_facet_orientation_info.flipped ⊻ neighbor_facet_orientation_info.flipped)
        shift = cell_facet_orientation_info.shift_index - (flip_interface ? neighbor_facet_orientation_info.shift_index : length(neighbor_children) - neighbor_facet_orientation_info.shift_index)
        new_neighbor_children = neighbor_children
        # URGERNT TODO: find a way to enable next line without allocations
        # new_neighbor_children = circshift((flip_interface ? reverse(neighbor_children) : neighbor_children), shift)
        new_neighbor_children = reverse(neighbor_children)
        # neighbor_local_seq = (neighbor_cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * (max(0, get_refinement_level(cell)))))) >>  ((Dim + 1) * (max(0, get_refinement_level(cell)) ))
        # @show
        if (get_refinement_level(neighbor_cell) > get_refinement_level(cell))
            neighbor_local_seq = (neighbor_cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * (get_refinement_level(neighbor_cell) - get_refinement_level(cell) -1)))) >>  ((Dim + 1) * (get_refinement_level(neighbor_cell) - get_refinement_level(cell) -1))
            # if get_refinement_level(cell) == 0
            #     neighbor_local_seq = (neighbor_cell.sequence & ((1 << (Dim + 1)) - 1))
            # end
            if findfirst(==(new_seq), cell_children) == findfirst(==(neighbor_local_seq), new_neighbor_children)
                if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor_idx]]
                    neighbors_count += length(neighbor_children)
                    continue
                else
                    neighbors_count += 1
                    continue
                end
            end
        elseif (get_refinement_level(neighbor_cell) < get_refinement_level(cell))
            _new_seq = (cell.sequence & (((1 << (Dim + 1)) - 1) << ((Dim + 1) * (max(get_refinement_level(neighbor_cell) - 1, 0))))) >>  ((Dim + 1) * (max(get_refinement_level(neighbor_cell)  - 1, 0)))
            # if findfirst(==(_new_seq), cell_children) == findfirst(==(neighbor_local_seq), new_neighbor_children)
                if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor_idx]]
                    neighbors_count += 1
                    continue
                else
                    neighbors_count += 1
                    continue
                end
            # end
        elseif (get_refinement_level(neighbor_cell) == get_refinement_level(cell))
            if refinement_cache.marked_for_refinement[refinement_cache.old_cell_to_new_cell_map[neighbor_idx]]
                neighbors_count += 1
                continue
            else
                neighbors_count += 1
                continue
            end
        end
    end
    return neighbors_count
end


function update_coarsened_dofs!(
    grid::KoppGrid{Dim},
    refinement_cache::KoppRefinementCache,
    dh::DofHandler,
    temp_celldofs::Vector{Int},
    temp_cell_dofs_offset::Vector{Int},
    temp_cell_to_subdofhandler::Vector{Int}) where Dim
    # dh.cell_dofs_offset .= 0
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        new == 0 && continue
        dh.cell_dofs_offset[new] = temp_cell_dofs_offset[old]
    end
    for (cell_idx, cell) in enumerate(grid.kopp_cells)
        refinement_cache.marked_for_coarsening[cell_idx] || continue
        for (other_cell_idx, other_cell) in enumerate(grid.kopp_cells)
            refinement_cache.old_cell_to_new_cell_map[other_cell_idx] == 0 && continue
            dh.cell_dofs_offset[other_cell_idx] <= temp_cell_dofs_offset[cell_idx] && continue
            dh.cell_dofs_offset[other_cell_idx] -= dh.subdofhandlers[1].ndofs_per_cell * (2^Dim)
            dh.ndofs -= dh.subdofhandlers[1].ndofs_per_cell * (2^Dim)
        end
    end
    for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        # copy old dofhander vectors
        new == 0 && continue
        dh.cell_dofs[dh.cell_dofs_offset[new] : dh.cell_dofs_offset[new] + dh.subdofhandlers[1].ndofs_per_cell - 1] .= @view temp_celldofs[temp_cell_dofs_offset[old] : temp_cell_dofs_offset[old] + dh.subdofhandlers[1].ndofs_per_cell - 1]
        dh.cell_to_subdofhandler[new] = temp_cell_to_subdofhandler[old]
        refinement_cache.marked_for_coarsening[new] || continue
        parent_offset = dh.cell_dofs_offset[new]
        parent_sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[new]]
        parent_dofs = @view dh.cell_dofs[parent_offset : parent_offset + parent_sdh.ndofs_per_cell - 1]
        for new_seq in 1 : 2^Dim
            child_idx = old + new_seq
            # For DG only
            offset = temp_cell_dofs_offset[child_idx]
            dofs = @view temp_celldofs[offset : offset + parent_sdh.ndofs_per_cell - 1]
            parent_dofs[new_seq] = dofs[new_seq]
        end
    end
end

@generated function get_child_weighted_dofs(::Type{T}) where {Dim, T <: Ferrite.RefHypercube{Dim}}
    ip = Lagrange{T, 1}()
    ntuple(seq -> (
    coords = Ferrite.reference_coordinates(ip);
    map!(node -> (node + coords[seq])/2, coords, coords);
    nshape_functions = length(coords);
    # TODO URGENT check if rows and cols are not switched
    SMatrix{nshape_functions, nshape_functions, Float64}(Ferrite.reference_shape_value(ip, coords[child_node], i) for i in 1:nshape_functions, child_node in 1:nshape_functions)),2^Dim)
end

function interpolate_solution!(refshape::Type{<:Ferrite.RefHypercube}, u::Vector{Float64}, dofs, parent_dofs, seq)
    u_parent = @view u[parent_dofs]
    weights = get_child_weighted_dofs(refshape)[seq]
    for i in 1:length(dofs)
        u[dofs[i]] = u_parent ⋅ (@view weights[:,i])
    end
end

@generated function get_ref_face_neighborhood(::Type{T}) where {Dim, T <: Ferrite.RefHypercube{Dim}}
    I = Tensors.diagm(Tensor{2,Dim}, ones(Dim))
    normals = Ferrite.weighted_normal.(Ref(I), T, 1:2*Dim)
    return ntuple(facet -> findfirst(x -> x ≈ -Ferrite.weighted_normal(I, T, facet), normals) , Val(Dim * 2))
end

# 1 alloc for hexahedron only IDK why hehe
@generated function get_ref_cell_neighborhood(::Type{T}) where {Dim, T <: Ferrite.RefHypercube{Dim}}
    return ntuple(cell ->
        ntuple(facet -> cell ∈ Ferrite.reference_facets(T)[facet] ? 0 :
            Ferrite.reference_facets(T)[facet][findfirst(facet_node -> any(facet_node .∈ filter(nodes ->  cell ∈ nodes, Ferrite.reference_edges(T))), Ferrite.reference_facets(T)[facet])], Val(Dim * 2)),
        Val(2^Dim))
end

function _process_seq(coords::NTuple{NNodes, NodeType}, seq::TInt) where {NodeType, TInt, NNodes}
    Dim = TInt(log2(NNodes)) # Only first order
    level::TInt = 0
    nbits::TInt = Dim + 1
    mask = (1 << nbits) - 1 # Does this need to be T too?
    maximum_level::TInt = sizeof(TInt)*8 ÷ nbits # Maybe use ÷ ?
    coord_new = MVector{NNodes, NodeType}(coords)
    coord_new_temp = MVector{NNodes, NodeType}(coords)
    @inbounds for level::TInt in maximum_level:-1:0 # There should be an easier way
        local_seq::TInt = (seq & (mask << (nbits * level)) ) >> (nbits * level )
        local_seq == 0 && continue
        @inbounds for (i, coord) in enumerate(coord_new_temp)
            coord_new[i] = NodeType((coord.x + coord_new_temp[local_seq].x)/2)
        end
        coord_new_temp .= coord_new
    end
    return Tuple(coord_new)
end

function _process_seq(Dim, coords::NTuple{NNodes, NodeType}, seq::TInt) where {NodeType, TInt, NNodes}
    level::TInt = 0
    nbits::TInt = Dim + 1
    mask = (1 << nbits) - 1 # Does this need to be T too?
    maximum_level::TInt = sizeof(TInt)*8 ÷ nbits # Maybe use ÷ ?
    coord_new = MVector{NNodes, NodeType}(coords)
    coord_new_temp = MVector{NNodes, NodeType}(coords)
    @inbounds for level::TInt in maximum_level:-1:0 # There should be an easier way
        local_seq::TInt = (seq & (mask << (nbits * level)) ) >> (nbits * level )
        local_seq == 0 && continue
        @inbounds for (i, coord) in enumerate(coord_new_temp)
            coord_new[i] = NodeType((coord.x + coord_new_temp[local_seq].x)/2)
        end
        coord_new_temp .= coord_new
    end
    return Tuple(coord_new)
end


function _transform_to_parent!(coords::AbstractVector, node_coords, seq::T) where {T}
    Dim = 2
    level::T = 0
    nbits::T = Dim + 1
    mask = (1 << nbits) - 1 # Does this need to be T too?
    maximum_level::T = sizeof(T)*8 ÷ nbits # Maybe use ÷ ?
    # coord_new = MVector{length(coords), eltype(coords)}(coords)
    for level::T in maximum_level:0 # There should be an easier way
        local_seq = seq & (mask << (nbits * level))
        local_seq == 0 && continue
        map!(node -> ((node + node_coords[local_seq])/2), coords, coords)
    end
    return nothing
end

function get_root_idx(grid::KoppGrid, i::Int)
    parent = grid.kopp_cells[i].parent
    # TODO: URGENT maxiter
    while parent > 0
        parent = grid.kopp_cells[parent].parent
    end
    return -parent
end

# Works only for 1st order geometric interpolation
# Has 0 allocations
function get_refined_coords(grid::KoppGrid{Dim}, i::Int) where Dim
    cell = grid.kopp_cells[i]
    root_idx = get_root_idx(grid, i)
    node_ids = Ferrite.get_node_ids(grid.base_grid.cells[root_idx])
    root_coords = ntuple(i -> ( grid.base_grid.nodes[node_ids[i]]), Val(2^(Dim)))
    return _process_seq(root_coords, cell.sequence)
end

# # 5 allocs
# function get_refinement_face_neighborhood(grid::KoppGrid{Dim}, topology::KoppTopology, i::Int) where {Dim}
#     cell = grid.kopp_cells[i]
#     seq = cell.sequence
#     refshape = getrefshape(getcelltype(grid.base_grid))
#     ref_neighborhood = get_ref_face_neighborhood(refshape)
#     mask = (1 << (Dim + 1)) - 1 # Does this need to be T too?
#     local_seq = seq & mask
#     inherit_from_parent = local_seq .∈ Ferrite.reference_facets(refshape)
#     parent_neighbors = topology.neighbors[cell.parent]
#     return ntuple(facet -> inherit_from_parent[facet] ? copy(parent_neighbors[facet]) : [ref_neighborhood[facet]] , Val(Dim * 2))
# end

# @generated function get_children_cells_from_facet(::Type{T}) where {Dim, T <: Ferrite.RefHypercube{Dim}}
#     return Ferrite.reference_facets(T)
# end

function Ferrite.getcoordinates!(dst::Vector{Vec{dim,T}}, grid::KoppGrid, cell::Int) where{dim, T}
    for i in 1:length(dst)
        dst[i] =  get_refined_coords(grid, cell)[i].x
    end
end

function Ferrite.getcoordinates!(dst::Vector{Vec{dim,T}}, grid::KoppGrid, cell::CellIndex) where{dim, T}
    return Ferrite.getcoordinates!(dst, grid, cell.idx)
end
