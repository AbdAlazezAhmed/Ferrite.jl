abstract type AbstractAdaptiveGrid{dim} <: AbstractGrid{dim} end
abstract type AbstractAdaptiveCell{refshape <: AbstractRefShape} <: AbstractCell{refshape} end
const refinable_interpolations = [
    Lagrange
]
function transform_refined(::Lagrange{RefQuadrilateral}, coord::Vec, i::Int)
    displacements = (Vec((-0.5,-0.5)),Vec((0.5,-0.5)),Vec((0.5,0.5)),Vec((-0.5,0.5)))
    return coord/2 + displacements[i]
end
for IP in (:(Lagrange{RefQuadrilateral, 1}),
    :(Lagrange{RefQuadrilateral, 2}),
    :(Lagrange{RefQuadrilateral, 3}),
    )
    @eval begin  
        @generated function get_shared_dofs(::$IP)
            ref = reference_coordinates($IP())
            j = length(ref) + 1
            res = Int[]
            for i in 1:4
                for x in transform_refined.(Ref($IP()), reference_coordinates($IP()), i)
                    temp = findfirst(u -> x ≈ u, unique(ref))
                    if temp === nothing
                        push!(res,  j)
                        j+=1
                    else
                        push!(res,  temp)
                    end
                    push!(ref, x)
                end
            end
            return res
        end
    end
end

struct KoppCell <: AbstractCell{AbstractRefShape{2}}
    secquence::Int # TODO: use two integers instead? alloc == bad
    parent::Int
    children::NTuple{4, Int}
    root::Int
    neighbors::NTuple{4, Pair{Int, Int}}
    isleaf::Bool
end

struct KoppRoot <: AbstractCell{AbstractRefShape{2}}
    children::Vector{KoppCell}
    refinement_index::ScalarWrapper{Int}
end

function get_refinement_level(root::KoppRoot, cell::KoppCell)
    parent = cell.parent
    i = 1
    while parent > 0
        i+=1
        parent = root.children[parent].parent
    end
    return i
end

mutable struct KoppGrid{G}
    base_grid::G
    kopp_roots::Vector{KoppRoot}
end

struct KoppCellCache{dim,X,G<:AbstractGrid,DH<:Union{AbstractDofHandler,Nothing}}
    flags::UpdateFlags
    grid::G
    # Pretty useless to store this since you have it already for the reinit! call, but
    # needed for the CellIterator(...) workflow since the user doesn't necessarily control
    # the loop order in the cell subset.
    cellid::ScalarWrapper{Int}
    nodes::Vector{Int}
    coords::Vector{X}
    dh::DH
    dofs::Vector{Int}
    dofs_bitmask::BitArray{dim}
end
# Quadrilateral
function generate_grid(C::Type{KoppCell}, nel::NTuple{2,Int}) where {T}
    grid =  generate_grid(Quadrilateral, nel)
    topology = ExclusiveTopology(grid)
    cells = KoppRoot[]
    for cell in 1:getncells(grid)
        neighborhood = tuple([isempty(face) ? -1 => 0 : face[1][1] => 0 for face in @view topology.edge_edge_neighbor[cell, :]]...)
        push!(cells, Ferrite.KoppRoot(KoppCell[
            KoppCell(
            0,
            -1,
            ntuple(j -> j, 4),
            cell,
            neighborhood,
            true
        )
        ],ScalarWrapper(0)))
    end
    return Ferrite.KoppGrid(grid, cells)
end
function get_refined_neighbor(grid::G, neighboring_cell::Pair{Int, Int}, current_face::FaceIndex) where G <: KoppGrid
    neighboring_cell[1] == -1 && return -1 => 0
    cell = grid.kopp_roots[neighboring_cell[1]]
    if cell.children |> isempty
        return neighboring_cell => 0
    elseif current_face ∈ (FaceIndex(1,1), FaceIndex(3,2))
        return neighboring_cell[1] =>  4
    elseif current_face ∈ (FaceIndex(1,4), FaceIndex(3,3))
        return neighboring_cell[1] =>  2
    elseif current_face ∈ (FaceIndex(2,1), FaceIndex(4,4))
        return neighboring_cell[1] =>  3
    elseif current_face ∈ (FaceIndex(2,2), FaceIndex(4,3))
        return neighboring_cell[1] =>  1
    end
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



function get_refined_neighbor(grid::G, parent_neighbot_idx::Pair{Int, Int}, current_face::FaceIndex) where G <: KoppGrid
    neighbor_root_idx = parent_neighbot_idx[1]
    neighbor_cell_idx = parent_neighbot_idx[2]
    neighbor_root_idx == -1 && return -1 => 0
    neighbor_root = grid.kopp_roots[neighbor_root_idx]
    neighbor_cell = neighbor_root.children[neighbor_cell_idx]
    @info neighbor_cell
    @info get_refinement_level(neighbor_root, neighbor_cell)

    if neighbor_cell.children == zeros(4)
        return neighboring_cell => 0
    elseif current_face ∈ (FaceIndex(1,1), FaceIndex(3,2))
        return neighboring_cell[1] =>  4
    elseif current_face ∈ (FaceIndex(1,4), FaceIndex(3,3))
        return neighboring_cell[1] =>  2
    elseif current_face ∈ (FaceIndex(2,1), FaceIndex(4,4))
        return neighboring_cell[1] =>  3
    elseif current_face ∈ (FaceIndex(2,2), FaceIndex(4,3))
        return neighboring_cell[1] =>  1
    end
end


function refine!(kopp_grid::KoppGrid, cell_idx::Pair{Int, Int})
    kopp_root = kopp_grid.kopp_roots[cell_idx[1]]
    parent = kopp_grid.kopp_roots[cell_idx[1]]
    idx = kopp_root.refinement_index[]
    refined_neighborhood_table = SMatrix{4,4}(
    0,2,4,0,
    0,0,3,1,
    2,0,0,4,
    1,3,0,0
    )
    for i in 1:4
        push!(kopp_root.children, KoppCell(
            i,
            cell_idx[2],
            ntuple(j -> 0, 4),
            # ntuple(j -> idx + j, 4),
            cell_idx[1],
            ntuple( j ->(cell_idx[1] => (idx + refined_neighborhood_table[i,j])) , 4),
            true
        ))
    end

    kopp_root.refinement_index[] += 4

end
