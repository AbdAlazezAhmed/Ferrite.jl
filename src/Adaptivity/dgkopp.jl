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
    refinement_order::Vector{Int}
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
        zeros(Int, 1)))
    end
    return Ferrite.DGKoppGrid(grid, cells)
end

function mark_for_refinement(grid::DGKoppGrid, cellset::Vector{Pair{Int, Int}})
    @assert all([all(root.refinement_order .== 0) for root in grid.kopp_roots]) "Grid must have no marked for refinement cells before marking for refinement"
    for (i, idx) in enumerate(cellset)
        root = grid.kopp_roots[idx[1]]
        root.refinement_order[idx[2]] = i
    end
end

function get_refined_neighbor(grid::G, neighboring_cell::Pair{Int, Int}, current_face::FaceIndex, n::Pair{Int, Int}) where G <: DGKoppGrid
    neighboring_cell[1] == -1 && return -1 => 0
    root = grid.kopp_roots[neighboring_cell[1]]
    cell = root.children[neighboring_cell[2] + 1]
    if cell.children == (0, 0, 0, 0)
        return neighboring_cell
    elseif current_face ∈ (FaceIndex(1,1), FaceIndex(3,2))
        cell_idx = cell.children[4]
    elseif current_face ∈ (FaceIndex(1,4), FaceIndex(3,3))
        cell_idx = cell.children[2]
    elseif current_face ∈ (FaceIndex(2,1), FaceIndex(4,4))
        cell_idx = cell.children[3]
    elseif current_face ∈ (FaceIndex(2,2), FaceIndex(4,3))
        cell_idx = cell.children[1]
    end
    cell = root.children[cell_idx + 1]
    neighboring_face = face_face_neighbor_local[current_face[2]]
    cell.neighbors = ntuple(i -> i == neighboring_face ? n : cell.neighbors[i], 4)
    return neighboring_cell[1] =>  cell_idx
end
refined_neighborhood_table = SMatrix{4,4}(
    0,2,4,0,
    0,0,3,1,
    2,0,0,4,
    1,3,0,0
    )
face_face_neighbor_local = SVector(3,4,1,2)
boundary_faces = SMatrix{4,2}(
    1,4,
    1,2,
    2,3,
    3,4
)


function get_inherited_neighbors(grid::DGKoppGrid, parent::DGKoppCell, face::FaceIndex)
    neighbors = Pair{Int, Int}[]
    for face_neighbors in parent.neighbors
        for neighbor_idx in face_neighbors
            neighbor_root = grid.kopp_roots[neighbor_idx[1]]
            neighbor = neighbor_root[neighbor_idx[2]]
            if !isempty(neighbor.secquence)
                if length(neighbor.secquence) < length(neighbor.secquence)
                    
                end
                
            end
        end
    end
end

function get_neighboring_cells(grid::DGKoppGrid, parent_cell_idx::Pair{Int, Int}, face::FaceIndex)
    parent_root = grid.kopp_roots[parent_cell_idx[1]]
    parent = parent_root[parent_cell_idx[2]]
    refined_neighborhood_table = SMatrix{4,4}(
    0,2,4,0,
    0,0,3,1,
    2,0,0,4,
    1,3,0,0
    )
    if  refined_neighborhood_table[face[2],face[1]] == 0
        neighbors = get_inherited_neighbors(grid, parent, face)
        correct_neighbors_neighbors(cell, neighbors)
    elseif  refined_neighborhood_table[face[2],face[1]] == 1
        parent_cell_idx
    else
        parent_cell_idx[1] => length(parent_root.children) + refined_neighborhood_table[face[2],face[1]]
    end
end

#Idea : ghost cells for derefinement?
function refine_cell(grid::DGKoppGrid, cell_idx::Pair{Int, Int})
    for i in 2:4
        neighborhood = ntuple(j -> get_neighboring_cells(grid, cell_idx, FaceIndex(i, j)), 4)
        push!(kopp_root.children, DGKoppCell(
            [i],
            neighborhood,
            (Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[], Pair{Int, Int}[]),
            Int[]
        ))
    end
end

