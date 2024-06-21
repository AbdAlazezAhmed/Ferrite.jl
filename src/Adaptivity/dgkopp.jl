struct NeighborIndex
    root::Int
    child::Int
    neighbor_face::Int
end

struct DGKoppCell
    secquence::Vector{Int}
    neighbors::NTuple{4, Vector{NeighborIndex}}
    interface_matrix_index::NTuple{4, Vector{Int}}
    dofs::Vector{Int} #predetermined size
end

struct DGKoppRoot{M <: AbstractMatrix}
    children::Vector{DGKoppCell}
    cell_matrices::Vector{M}
    refinement_order::Vector{Int} #TODO: replace with Vecotr{Bool}
    children_updated_indices::Vector{Int}
end

struct DGKoppGrid{G, M <: AbstractMatrix}
    base_grid::G
    kopp_roots::Vector{DGKoppRoot{M}}
    interface_matrices::Vector{M}
    root_face_orientation_info::Vector{NTuple{4, OrientationInfo}}
end

function generate_grid(C::Type{DGKoppCell}, nel::NTuple{2,Int}) where {T}
    grid =  generate_grid(Quadrilateral, nel)
    return DGKoppGrid(grid)
end

function DGKoppGrid(grid::G) where {G <: AbstractGrid}
    topology = ExclusiveTopology(grid)
    cells = DGKoppRoot{Matrix{Float64}}[]
    for cell in 1:getncells(grid)
        neighborhood = tuple([isempty(face) ? NeighborIndex[] : [NeighborIndex(face[1][1], 1, face[1][2])] for face in @view topology.edge_edge_neighbor[cell, :]]...)
        push!(cells, Ferrite.DGKoppRoot(DGKoppCell[
            DGKoppCell(
            Int[],
            neighborhood,
            (Int[], Int[], Int[], Int[]),
            Int[]        
            )
        ],
        Matrix{Float64}[],
        zeros(Int, 1),
        [1]))
    end
    interface_matrices = Matrix{Float64}[]
    resize!(interface_matrices, length(facetskeleton(topology, grid)))
    for (i, edge) in enumerate(facetskeleton(topology, grid))
        push!(cells[edge[1]].children[1].interface_matrix_index[edge[2]], i)
        neighbors_vec = cells[edge[1]].children[1].neighbors[edge[2]]
        isempty(neighbors_vec) && continue
        neighbor = first(neighbors_vec)
        push!(cells[neighbor.root].children[neighbor.child].interface_matrix_index[neighbor.neighbor_face], i)
    end
    orientation = [ntuple(i -> OrientationInfo(facets(cell)[i]), 4) for cell in grid.cells]
    return Ferrite.DGKoppGrid(grid, cells, interface_matrices, orientation)
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
        resize!(root.cell_matrices, l)
        root.children .= similar(temp, l)
        for (old, new) in enumerate(root.children_updated_indices)
            root.children[new] = deepcopy(temp[old])
            for j in 1:4
                for (neighbor_idx, neighbor) in enumerate(root.children[new].neighbors[j])
                    root.children[new].neighbors[j][neighbor_idx] = NeighborIndex(neighbor.root, grid.kopp_roots[neighbor.root].children_updated_indices[neighbor.child], neighbor.neighbor_face)
                end
            end
        end
    end
end

function get_inherited_neighbors(grid::DGKoppGrid, parent_cell_idx::Pair{Int, Int}, face::FaceIndex)
    parent = grid.kopp_roots[parent_cell_idx[1]].children[parent_cell_idx[2]]
    neighbors = NeighborIndex[]
    local_seq = face[1]
    local_face = face[2]
    refined_neighborhood_table = (
        SMatrix{4,4}([4 0 0 2
        3 1 0 0
        0 4 2 0
        0 0 1 3]),

        SMatrix{4,4}([3 0 0 1
        2 4 0 0
        0 3 1 0
        0 0 4 2]),

        SMatrix{4,4}([2 0 0 4
        1 3 0 0
        0 2 4 0
        0 0 3 1]),

        SMatrix{4,4}([1 0 0 3
        4 2 0 0
        0 1 3 0
        0 0 2 4])

        )

        
    _flipped = (3,4,1,2)
    __flipped = [
        2 0 0 4
        1 3 0 0
        0 2 4 0
        0 0 3 1
    ]
    neighbor_child_face_map = SMatrix{4,2}(
        1,2,3,4,
        2,3,4,1
    )
    for neighbor_idx in parent.neighbors[local_face]
        neighbor_root = grid.kopp_roots[neighbor_idx.root]
        neighbor = neighbor_root.children[neighbor_idx.child]
        idx = 1
        flipped = false
        if neighbor_idx.root != parent_cell_idx[1]
            idx = _flipped[face[2]] - neighbor_idx.neighbor_face
            idx < 0 && (idx += 4)
            idx += 1
            flipped = grid.root_face_orientation_info[neighbor_idx.root][neighbor_idx.neighbor_face].flipped  == grid.root_face_orientation_info[parent_cell_idx[1]][face[2]].flipped
        end
        if neighbor_root.refinement_order[neighbor_idx.child] == 0
            if !isempty(neighbor.secquence)
                if length(neighbor.secquence) > length(parent.secquence)
                    neighbor_seq = neighbor.secquence[length(parent.secquence) + 1]
                    if neighbor_seq == refined_neighborhood_table[idx][local_seq, local_face]
                        push!(neighbors, neighbor_idx)
                    end
                else
                    push!(neighbors, neighbor_idx) #Should be the only neighbor?
                end
            else
                push!(neighbors, neighbor_idx) #Should be the only neighbor?
            end
        else # The parent neighbor is marked for refinement 
            if length(neighbor.secquence) > length(parent.secquence)
                if neighbor.secquence[length(parent.secquence) + 1] == refined_neighborhood_table[idx][local_seq, local_face]
                    _irange = flipped ? (2:-1:1) : 1:2
                    for i in _irange
                        push!(neighbors, NeighborIndex(neighbor_idx.root, neighbor_idx.child + neighbor_child_face_map[neighbor_idx.neighbor_face, i] - 1, neighbor_idx.neighbor_face))
                    end
                end
            elseif length(neighbor.secquence) == length(parent.secquence)
                push!(neighbors, NeighborIndex(neighbor_idx.root, neighbor_idx.child + (flipped ? __flipped[refined_neighborhood_table[idx][local_seq, local_face], neighbor_idx.neighbor_face] : refined_neighborhood_table[idx][local_seq, local_face]) - 1, neighbor_idx.neighbor_face))
            end
        end
    end
    return neighbors
end

function correct_neighbors_neighbors!(grid::DGKoppGrid, cell::Pair{Int, Int}, neighbors::Vector{NeighborIndex}, face::FaceIndex)
    face_face_neighbor_local = (3,4,1,2)
    for neighbor_idx in neighbors
        neighbor_root = grid.kopp_roots[neighbor_idx.root]
        isassigned(neighbor_root.refinement_order, neighbor_idx.child) || continue
        neighbor_root.refinement_order[neighbor_idx.child] != 0 && continue
        neighbor_cell = neighbor_root.children[neighbor_idx.child]
        neighbor_face = neighbor_idx.neighbor_face
        filter!(x -> x != NeighborIndex(cell[1], cell[2] - face[1] + 1, face[2]), neighbor_cell.neighbors[neighbor_face])
        push!(neighbor_cell.neighbors[neighbor_face], NeighborIndex(cell[1], cell[2], face[2]))
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
    local_face_face = (3,4,1,2)
    if  refined_neighborhood_table[face[2],face[1]] == 0
        neighbors = get_inherited_neighbors(grid, parent_cell_idx, face)
        return neighbors
    elseif  refined_neighborhood_table[face[2],face[1]] == 1
        return [NeighborIndex(parent_cell_idx[1], parent_cell_idx[2], local_face_face[face[2]])]
    else
        return [NeighborIndex(parent_cell_idx[1], parent_cell_idx[2] + refined_neighborhood_table[face[2],face[1]] - 1, local_face_face[face[2]])]
    end
end

function refine!(grid::DGKoppGrid)
    neighborhoods = NTuple{4, Vector{NeighborIndex}}[]
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
        root = grid.kopp_roots[root_idx]
        old_cell_idx = 1
        for new_cell_idx in root.children_updated_indices
            isassigned(root.children, new_cell_idx) || continue
            if root.refinement_order[old_cell_idx] != 0
                for _i in 1:4
                    neighborhood = neighborhoods[current_cell_idx]
                    for k in eachindex(neighborhood)
                        correct_neighbors_neighbors!(grid, root_idx => new_cell_idx + _i - 1, neighborhood[k], FaceIndex(_i, k))
                    end
                    root.children[new_cell_idx + _i - 1] = DGKoppCell(
                        [_i],
                        neighborhood,
                        (Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[]),
                        Int[]
                    )
                    current_cell_idx += 1
                end
            end
            old_cell_idx += 1
        end
    end
end

# Works only for 1st order geometric interpolation
function transform_refined(coord::NTuple{4, <:Vec}, i::Int)
    return ntuple(j -> (coord[i]+ coord[j])/2, 4)
end

function get_coords(grid::DGKoppGrid, cell_idx::Pair{Int, Int})
    root = grid.kopp_roots[cell_idx[1]]
    cell = root.children[cell_idx[2]]
    nodes_idx = grid.base_grid.cells[cell_idx[1]].nodes
    coords = ntuple(i -> grid.base_grid.nodes[nodes_idx[i]].x, 4)
    for k in length(cell.secquence):-1:1
        coords = transform_refined(coords, cell.secquence[k])
    end
    return coords
end

# 2 allocs for J?
function assemble_element_matrix!(ip::DiscontinuousLagrange, qr::QuadratureRule, grid::DGKoppGrid, cell::Pair{Int, Int})
    root = grid.kopp_roots[cell[1]]
    n_basefuncs = getnbasefunctions(ip)
    if !isassigned(root.cell_matrices, cell[2])
        root.cell_matrices[cell[2]] = zeros(n_basefuncs, n_basefuncs)
    else
        fill!(root.cell_matrices[cell[2]], 0)
    end
    Ke = root.cell_matrices[cell[2]]
    geo_mapping = GeometryMapping{1}(Float64, ip, qr)
    x = get_coords(grid, cell)
    for q_point in 1:getnquadpoints(qr)
        ξ = qr.points[q_point]
        J = calculate_mapping(geo_mapping, q_point, x).J
        w = qr.weights[q_point]
        dΩ = calculate_detJ(J) * w
        for i in 1:n_basefuncs
            δu  = shape_value(ip, ξ, i)
            for j in 1:n_basefuncs
                u = shape_value(ip, ξ, i)
                Ke[i, j] += (δu ⋅ u) * dΩ
            end
        end
    end
    return Ke
end


function assemble_interface_matrix!(ip::DiscontinuousLagrange, qr::FacetQuadratureRule, grid::DGKoppGrid, cell::Pair{Int, Int})
    root = grid.kopp_roots[cell[1]]
    n_basefuncs = getnbasefunctions(ip)
    if !isassigned(root.cell_matrices, cell[2])
        root.cell_matrices[cell[2]] = zeros(n_basefuncs, n_basefuncs)
    else
        fill!(root.cell_matrices[cell[2]], 0)
    end
    Ke = root.cell_matrices[cell[2]]
    geo_mapping = GeometryMapping{1}(Float64, ip, qr)
    x = get_coords(grid, cell)
    for q_point in 1:getnquadpoints(qr)
        ξ = qr.points[q_point]
        J = calculate_mapping(geo_mapping, q_point, x).J
        w = qr.weights[q_point]
        dΩ = calculate_detJ(J) * w
        for i in 1:n_basefuncs
            δu  = shape_value(ip, ξ, i)
            for j in 1:n_basefuncs
                u = shape_value(ip, ξ, i)
                Ke[i, j] += (δu ⋅ u) * dΩ
            end
        end
    end
    return Ke
end