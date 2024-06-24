struct NeighborIndex
    root::Int
    child::Int
    neighbor_face::Int
end

struct DGKoppCell
    secquence::Vector{Int}
    neighbors::NTuple{4, Vector{NeighborIndex}}
    interface_matrix_index::NTuple{4, Vector{Int}}
end

struct DGKoppRoot{M <: AbstractMatrix}
    children::Vector{DGKoppCell}
    cell_matrices::Vector{M}
    marked_for_refinement::BitVector
    dofs::Vector{Int}
    children_updated_indices::Vector{Int} #TODO: remove?
end

struct DGKoppGrid{G, M <: AbstractMatrix}
    base_grid::G
    kopp_roots::Vector{DGKoppRoot{M}}
    interface_matrices::Vector{M}
    root_face_orientation_info::Vector{NTuple{4, OrientationInfo}}
    interfaces_updated_indices::Vector{Int} #TODO: remove?
    interfaces_recompute::BitVector
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
            (Int[], Int[], Int[], Int[])     
            )
        ],
        Matrix{Float64}[],
        falses(1),
        zeros(4),
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
    return Ferrite.DGKoppGrid(grid, cells, interface_matrices, orientation, collect(1:length(interface_matrices)), falses(length(interface_matrices)))
end

function mark_for_refinement(grid::DGKoppGrid, cellset::Set{Pair{Int, Int}})
    @assert all([all(.!root.marked_for_refinement) for root in grid.kopp_roots]) "Grid must have no marked for refinement cells before marking for refinement"
    refined_interfaces = Int[]
    for (i, root) in enumerate(grid.kopp_roots)
        k = 0
        k_i = 0
        for (j, cell) in enumerate(root.children)
            if (i => j) ∈ cellset
                k += 1
                k_i += 1
                root.marked_for_refinement[j] = true
                root.children_updated_indices[j+1:end] .+= 3
                last_iterated_interface = 0
                for face in cell.interface_matrix_index
                    for interface in face
                        interface ∈ refined_interfaces && continue
                        grid.interfaces_recompute[interface:interface + 1] .= true
                        grid.interfaces_updated_indices[interface+1:end] .+=1
                        last_iterated_interface = interface
                        push!(refined_interfaces, interface)
                    end
                end
                if last_iterated_interface != 0
                    grid.interfaces_updated_indices[last_iterated_interface+1:end] .+=4
                end
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
    resize!(grid.interface_matrices, length(grid.interface_matrices) + length(refined_interfaces) + 4*length(cellset))

    return nothing
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
        if neighbor_root.marked_for_refinement[neighbor_idx.child] == false
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
        isassigned(neighbor_root.marked_for_refinement, neighbor_idx.child) || continue
        neighbor_root.marked_for_refinement[neighbor_idx.child] && continue
        neighbor_cell = neighbor_root.children[neighbor_idx.child]
        neighbor_face = neighbor_idx.neighbor_face
        @info "cells before delete" neighbor_cell.neighbors[neighbor_face]
        filter!(x -> x != NeighborIndex(cell[1], cell[2] - face[1] + 1, face[2]), neighbor_cell.neighbors[neighbor_face])
        @info "cells after delete" neighbor_cell.neighbors[neighbor_face]
        push!(neighbor_cell.neighbors[neighbor_face], NeighborIndex(cell[1], cell[2], face[2]))
        @info "cells after add" neighbor_cell.neighbors[neighbor_face]
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
        for (old_cell_idx, new_cell_idx) in enumerate(root.children_updated_indices)
            isassigned(root.children, new_cell_idx) || continue
            if root.marked_for_refinement[old_cell_idx]
                for _i in 4:-1:1 # hotfix
                    neighborhood = ntuple(k -> get_neighboring_cells(grid, root_idx=>new_cell_idx, FaceIndex(_i, k)), 4)
                    push!(neighborhoods, neighborhood)
                end
            else
                neighborhood = root.children[new_cell_idx].neighbors
                push!(neighborhoods, neighborhood)
            end

        end
    end
    @info neighborhoods
    current_cell_idx = 1
    for root_idx in eachindex(grid.kopp_roots)
        root = grid.kopp_roots[root_idx]
        old_cell_idx = 1
        for new_cell_idx in root.children_updated_indices
            isassigned(root.children, new_cell_idx) || continue
            if root.marked_for_refinement[old_cell_idx]
                for _i in 4:-1:1 # hotfix
                    neighborhood = neighborhoods[current_cell_idx]
                    for k in eachindex(neighborhood)
                        correct_neighbors_neighbors!(grid, root_idx => new_cell_idx + _i - 1, neighborhood[k], FaceIndex(_i, k))
                    end
                    parent = root.children[new_cell_idx]
                    root.children[new_cell_idx + _i - 1] = DGKoppCell(
                        [parent.secquence; _i],
                        neighborhood,
                        (Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[]),
                    )
                    current_cell_idx += 1
                end
            else
                current_cell_idx += 1
            end
            old_cell_idx += 1
        end
        resize!(grid.kopp_roots[root_idx].marked_for_refinement, length(grid.kopp_roots[root_idx].children))
        grid.kopp_roots[root_idx].marked_for_refinement .= falses(length(grid.kopp_roots[root_idx].marked_for_refinement))
        resize!(grid.kopp_roots[root_idx].children_updated_indices, length(grid.kopp_roots[root_idx].children))
        grid.kopp_roots[root_idx].children_updated_indices .= collect(1:length(grid.kopp_roots[root_idx].children))
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


function assemble_interface_matrix!(ip::DiscontinuousLagrange, qr::FacetQuadratureRule, grid::DGKoppGrid, cell_idx::Pair{Int, Int}, face::Int)
    root = grid.kopp_roots[cell_idx[1]]
    cell = root.children[cell_idx[2]]
    n_basefuncs = getnbasefunctions(ip)
    for (i, i_matrix) in enumerate(cell.interface_matrix_index[face])
        if cell.neighbors[face] |> isempty
            return nothing
        else
            neighbor_idx = cell.neighbors[face][i]
            neighbor_root = grid.kopp_roots[neighbor_idx.root]
            neighbor_cell = neighbor_root.children[neighbor_idx.child]
            cell_a = length(neighbor_cell.secquence) > length(cell.secquence) ? cell : neighbor_cell
            cell_a_idx = length(neighbor_cell.secquence) > length(cell.secquence) ? cell_idx : (neighbor_idx.root => neighbor_idx.child)
            cell_b = length(neighbor_cell.secquence) > length(cell.secquence) ? neighbor_cell : cell
            cell_b_idx = length(neighbor_cell.secquence) > length(cell.secquence) ? (neighbor_idx.root => neighbor_idx.child) : cell_idx
            face_a = length(neighbor_cell.secquence) > length(cell.secquence) ? face : neighbor_idx.neighbor_face
            face_b = length(neighbor_cell.secquence) > length(cell.secquence) ? neighbor_idx.neighbor_face : face

            OI_a = grid.root_face_orientation_info[cell_a_idx[1]][face_a]
            OI_b = grid.root_face_orientation_info[cell_b_idx[1]][face_b]
            flipped = OI_a.flipped != OI_b.flipped
            shift_index = OI_b.shift_index - OI_a.shift_index
            IOI = InterfaceOrientationInfo{RefQuadrilateral, RefQuadrilateral}(flipped, shift_index, OI_b.shift_index, face_a, face_b)

            if !isassigned(grid.interface_matrices, i_matrix)
                grid.interface_matrices[i_matrix] = zeros(2*n_basefuncs, 2*n_basefuncs)
            else
                fill!(grid.interface_matrices[i_matrix], 0.0)
            end
            Ki = grid.interface_matrices[i_matrix]
            geo_mapping_a = GeometryMapping{1}(Float64, ip, qr.face_rules[face_a])
            geo_mapping_b = GeometryMapping{1}(Float64, ip, qr.face_rules[face_b])
            x_a = get_coords(grid, cell_idx)
            x_b = get_coords(grid, neighbor_idx.root => neighbor_idx.child)
            transform_interface_points!
            for q_point in 1:getnquadpoints(qr.face_rules[face_a])
                ξₐ = qr.face_rules[face_a].points[q_point]
                face_point = element_to_facet_transformation(ξₐ, RefQuadrilateral, face_a)
                flipped && (face_point *= -1)
                ξᵦ = facet_to_element_transformation(face_point, RefQuadrilateral, face_b)
                Jₐ = calculate_mapping(geo_mapping_a, q_point, x_a).J
                Jᵦ = calculate_mapping(geo_mapping_b, q_point, x_b).J
                w = qr.face_rules[face_a].weights[q_point]
                Jₐ_inv = calculate_Jinv(Jₐ)
                Jᵦ_inv = calculate_Jinv(Jᵦ)
                weight_norm_a = weighted_normal(Jₐ, RefQuadrilateral, face_a)
                weight_norm_b = weighted_normal(Jᵦ, RefQuadrilateral, face_b)
                dΩₐ = norm(weight_norm_a) * w
                dΩᵦ = norm(weight_norm_b) * w
                dΓ = dΩₐ
                normal =  weight_norm_a / norm(weight_norm_a)
                for i in 1:2*n_basefuncs
                    ∇δu_here  = i > n_basefuncs ? shape_gradient(ip, i > n_basefuncs ? ξᵦ : ξₐ, i > n_basefuncs ? i - n_basefuncs  : i) : Vec(0.,0.)
                    ∇δu_there  = i > n_basefuncs ? Vec(0.,0.) : shape_gradient(ip, i > n_basefuncs ? ξᵦ : ξₐ, i > n_basefuncs ? i - n_basefuncs  : i)
                    δu_here  = i > n_basefuncs ? shape_value(ip, i > n_basefuncs ? ξᵦ : ξₐ, i > n_basefuncs ? i - n_basefuncs  : i) : 0.0
                    δu_there  = i > n_basefuncs ? 0.0 : shape_value(ip, i > n_basefuncs ? ξᵦ : ξₐ, i > n_basefuncs ? i - n_basefuncs  : i)
                    test_jump = -(δu_here - δu_there) * normal
                    test_grad_avg = (dothelper(∇δu_here, Jₐ_inv)+ dothelper(∇δu_there, Jᵦ_inv))/2
                    for j in 1:2*n_basefuncs 
                        ∇u_here  = j > n_basefuncs ? shape_gradient(ip, j > n_basefuncs ? ξᵦ : ξₐ, j > n_basefuncs ? j - n_basefuncs  : j) : Vec(0.,0.)
                        ∇u_there  = j > n_basefuncs ? Vec(0.,0.) : shape_gradient(ip, j > n_basefuncs ? ξᵦ : ξₐ, j > n_basefuncs ? j - n_basefuncs  : j)
                        u_here  = j > n_basefuncs ? shape_value(ip, j > n_basefuncs ? ξᵦ : ξₐ, j > n_basefuncs ? j - n_basefuncs  : j) : 0.0
                        u_there  = j > n_basefuncs ? 0.0 : shape_value(ip, j > n_basefuncs ? ξᵦ : ξₐ, j > n_basefuncs ? j - n_basefuncs  : j)
                        trial_jump = -(u_here - u_there) * normal
                        trial_grad_avg = (dothelper(∇u_here, Jₐ_inv)+ dothelper(∇u_there, Jᵦ_inv))/2
                        Ki[i, j] += -(test_jump ⋅ trial_grad_avg + test_grad_avg ⋅ trial_jump) * dΓ + 2 * (test_jump ⋅ trial_jump) * dΓ
                    end
                end
            end
        end
    end
end