struct DGKoppCell
    secquence::Vector{Int}
    neighbors::NTuple{4, Vector{Pair{Int, Int}}}
    interface_matrix_index::NTuple{4, Vector{Pair{Int, Int}}}
    dofs::Vector{Int} #predetermined size
end

struct DGKoppRoot{M <: AbstractMatrix}
    children::Vector{DGKoppCell}
    cell_matrices::Vector{M}
    interface_matrices::Vector{M}
    refinement_order::Vector{Int} #TODO: replace with Vecotr{Bool}
    children_updated_indices::Vector{Int}
end

struct DGKoppGrid{G}
    base_grid::G
    kopp_roots::Vector{DGKoppRoot}
end

function generate_grid(C::Type{DGKoppCell}, nel::NTuple{2,Int}) where {T}
    grid =  generate_grid(Quadrilateral, nel)
    topology = ExclusiveTopology(grid)
    cells = DGKoppRoot[]
    for cell in 1:getncells(grid)
        neighborhood = tuple([isempty(face) ? Pair{Int, Int}[] : [face[1][1] => 1] for face in @view topology.edge_edge_neighbor[cell, :]]...)
        push!(cells, Ferrite.DGKoppRoot(DGKoppCell[
            DGKoppCell(
            Int[],
            neighborhood,
            (Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[]),
            Int[]        
            )
        ],
        Matrix{Float64}[],
        Matrix{Float64}[],
        zeros(Int, 1),
        [1]))
    end
    return Ferrite.DGKoppGrid(grid, cells)
end

function mark_for_refinement(grid::DGKoppGrid, cellset::Set{Pair{Int, Int}})
    @assert all([all(root.refinement_order .== 0) for root in grid.kopp_roots]) "Grid must have no marked for refinement cells before marking for refinement"
    for (i, root) in enumerate(grid.kopp_roots)
        k = 0
        for (j, cell) in enumerate(root.children)
            if (i => j) ∈ cellset
                k+=1
                root.refinement_order[j] = k
                root.children_updated_indices[j+1:end] .+= 3
            end
        end
        temp = deepcopy(root.children)
        l = length(root.children)+3*k
        resize!(root.children, l)
        root.children .= similar(temp, l)
        for (old, new) in enumerate(root.children_updated_indices)
            root.children[new] = deepcopy(temp[old])
            for j in 1:4
                for (neighbor_idx, neighbor) in enumerate(root.children[new].neighbors[j])
                    root.children[new].neighbors[j][neighbor_idx] = neighbor[1] => grid.kopp_roots[neighbor[1]].children_updated_indices[neighbor[2]]
                end
            end
        end
    end
end

function get_inherited_neighbors(grid::DGKoppGrid, parent_cell_idx::Pair{Int, Int}, face::FaceIndex)
    parent = grid.kopp_roots[parent_cell_idx[1]].children[parent_cell_idx[2]]
    neighbors = Pair{Int, Int}[]
    local_seq = face[1]
    local_face = face[2]
    refined_neighborhood_table = SMatrix{4,4}(
        4,0,0,2,
        3,1,0,0,
        0,4,2,0,
        0,0,1,3
        )
    neighbor_child_face_map = SMatrix{4,2}(
        4,3,4,1,
        2,1,2,3
    )
    for neighbor_idx in parent.neighbors[local_face]
        neighbor_root = grid.kopp_roots[neighbor_idx[1]]
        neighbor = neighbor_root.children[neighbor_idx[2]]
        if neighbor_root.refinement_order[neighbor_idx[2]] == 0
            if !isempty(neighbor.secquence)
                if length(neighbor.secquence) > length(parent.secquence)
                    if neighbor.secquence[length(parent.secquence) + 1] == refined_neighborhood_table[local_face, local_seq]
                        push!(neighbors, neighbor_idx)
                    end
                else
                    push!(neighbors, neighbor_idx) #Should be the only neighbor?
                end
            else
                push!(neighbors, neighbor_idx) #Should be the only neighbor?
            end
            correct_neighbors_neighbors!(grid, parent_cell_idx[1] => parent_cell_idx[2] + face[1] - 1, neighbors, face)
        else # The parent neighbor is marked for refinement 
            if length(neighbor.secquence) > length(parent.secquence)
                if neighbor.secquence[length(parent.secquence) + 1] == refined_neighborhood_table[local_face, local_seq]
                    for i in 1:2
                        push!(neighbors, neighbor_idx[1] => neighbor_idx[2] + neighbor_child_face_map[local_seq, i] - 1)
                    end
                end
            elseif length(neighbor.secquence) == length(parent.secquence)
                push!(neighbors, neighbor_idx[1] => neighbor_idx[2] + refined_neighborhood_table[face[2], local_seq] - 1)
            end
        end
    end
    return neighbors
end

function correct_neighbors_neighbors!(grid::DGKoppGrid, cell::Pair{Int, Int}, neighbors::Vector{Pair{Int, Int}}, face::FaceIndex)
    face_face_neighbor_local = (3,4,1,2)
    for neighbor_idx in neighbors
        neighbor_root = grid.kopp_roots[neighbor_idx[1]]
        neighbor_root.refinement_order[neighbor_idx[2]] != 0 && continue
        neighbor_cell = neighbor_root.children[neighbor_idx[2]]
        neighbor_face = face_face_neighbor_local[face[2]]
        filter!(x -> x != (cell[1] => cell[2] - face[1] + 1), neighbor_cell.neighbors[neighbor_face])
        push!(neighbor_cell.neighbors[neighbor_face], cell)
    end
end

function get_neighboring_cells(grid::DGKoppGrid, parent_cell_idx::Pair{Int, Int}, face::FaceIndex)
    parent_root = grid.kopp_roots[parent_cell_idx[1]]
    parent = parent_root.children[parent_cell_idx[2]]
    refined_neighborhood_table = SMatrix{4,4}(
    0,2,4,0,
    0,0,3,1,
    2,0,0,4,
    1,3,0,0
    )
    if  refined_neighborhood_table[face[2],face[1]] == 0
        neighbors = get_inherited_neighbors(grid, parent_cell_idx, face)
        return neighbors
    elseif  refined_neighborhood_table[face[2],face[1]] == 1
        return [parent_cell_idx]
    else
        return [parent_cell_idx[1] => parent_cell_idx[2] + refined_neighborhood_table[face[2],face[1]] - 1]
    end
end

function refine!(grid::DGKoppGrid)
    neighborhoods = NTuple{4, Vector{Pair{Int, Int}}}[]
    for root_idx in eachindex(grid.kopp_roots)
        root = grid.kopp_roots[root_idx]
        for new_cell_idx in root.children_updated_indices
            isassigned(root.children, new_cell_idx) || continue
            for _i in 1:4
                neighborhood = ntuple(k -> get_neighboring_cells(grid, root_idx=>new_cell_idx, FaceIndex(_i, k)), 4)
                push!(neighborhoods, neighborhood)
            end
        end
    end
    current_cell_idx = 1
    for root_idx in eachindex(grid.kopp_roots)
        # Threads.@threads for root_idx in eachindex(grid.kopp_roots)
        root = grid.kopp_roots[root_idx]
        old_cell_idx = 1
        for new_cell_idx in root.children_updated_indices
            isassigned(root.children, new_cell_idx) || continue
            if root.refinement_order[old_cell_idx] != 0
                for _i in 1:4
                    neighborhood = neighborhoods[current_cell_idx]
                    root.children[new_cell_idx + _i - 1] = DGKoppCell(
                        [_i],
                        neighborhood,
                        (Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[]),
                        Int[]
                    )
                    current_cell_idx += 1
                end
                #TODO: modify the parent cell to be cell 1 in refinement
            end
            old_cell_idx += 1
        end
    end
end
