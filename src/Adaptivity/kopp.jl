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

# 
# struct KoppCell <: AbstractCell{AbstractRefShape{2}, GArray <: AbstractArray}
    # parent_secquence::GArray{Pair{Int, Int}}
    # root::Int
    # neighbors::GArray{4, Int}
    # isleaf::Bool
# end
# 

struct KoppCell <: AbstractCell{AbstractRefShape{2}}
    parent_secquence::Vector{Int} # TODO: use two integers instead? alloc == bad
    children::MVector{4, Int}
    root::Int
    neighbors::MVector{4, Pair{Int, Int}}
    isleaf::Bool
end

struct KoppRoot <: AbstractCell{AbstractRefShape{2}}
    children::Vector{KoppCell}
    refinement_index::ScalarWrapper{Int}
    neighbors::MVector{4, Int}
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
        neighborhood = MVector{4}([isempty(face) ? -1 : face[1][1] for face in topology.edge_edge_neighbor[cell, :]])
        push!(cells, Ferrite.KoppRoot(KoppCell[],ScalarWrapper(0), neighborhood))
    end
    return Ferrite.KoppGrid(grid, cells)
end
function get_refined_neighbor(grid::KoppGrid, neighboring_cell::Int, current_face::FaceIndex)
    neighboring_cell == -1 && return -1 => 0
    cell = grid.kopp_roots[neighboring_cell]
    if cell.children |> isempty
        return neighboring_cell => 0
    elseif current_face ∈ (FaceIndex(1,1), FaceIndex(3,2))
        return neighboring_cell =>  4
    elseif current_face ∈ (FaceIndex(1,4), FaceIndex(3,3))
        return neighboring_cell =>  2
    elseif current_face ∈ (FaceIndex(2,1), FaceIndex(4,4))
        return neighboring_cell =>  3
    elseif current_face ∈ (FaceIndex(2,2), FaceIndex(4,3))
        return neighboring_cell =>  1
    end
end
refined_neighborhood_table = (
    (0,2,4,0),
    (0,0,3,1),
    (2,0,0,4),
    (1,3,0,0)
    )
face_face_neighbor_local = (3,4,1,2)
boundary_faces = (
    (1,4),
    (1,2),
    (2,3),
    (3,4)
)
function refine_root!(kopp_grid::KoppGrid, cell_idx::Int)
    kopp_root = kopp_grid.kopp_roots[cell_idx]
    parent = kopp_grid.kopp_roots[cell_idx]
    append!(kopp_root.children,
        ntuple(i -> KoppCell(
            [i],
            MVector{4}(1,2,3,4),
            cell_idx,
            MVector{4}(ntuple(j -> refined_neighborhood_table[i][j] == 0 ? Ferrite.get_refined_neighbor(kopp_grid,parent.neighbors[j],FaceIndex(i,j)) : (cell_idx => (parent.refinement_index[] + refined_neighborhood_table[i][j])), 4)),
            true
        ), 4)
    )
    parent.refinement_index[] += 4
    # struct KoppCell <: AbstractCell{AbstractRefShape{2}}
    #     parent_secquence::Vector{Pair{Int, Int}}
    #     root::Int
    #     neighbors::MVector{4, Int}
    #     isleaf::Bool
    # end
end
function refine!(kopp_grid::KoppGrid, cell::Int)
    nkoppcells = length(kopp_grid.kopp_cells)
    parent = kopp_grid.kopp_cells[cell]
    refinement_level = parent.refinement_level + 1
    for i in 1:4
        parent_secquence = Int[]
        children = Int[0,0,0,0]
        rnt = refined_neighborhood_table[i]
        neighbors = [rnt[j] == 0 ? Ferrite.get_refined_neighbor(kopp_grid,parent.neighbors[j],FaceIndex(i,j)) : nkoppcells + rnt[j] for j in 1:4]
        kc = KoppCell(parent_secquence, refinement_level, cell, parent.root, neighbors, children)
        push!(kopp_grid.kopp_cells, kc)

    end
    for i in 1:4
        for (jj, j) in enumerate(kopp_grid.kopp_cells[nkoppcells + i].neighbors)
            j == -1  && continue
            neighbor = kopp_grid.kopp_cells[j]
            neighbor.refinement_level != refinement_level && continue
            neighbor.neighbors[face_face_neighbor_local[jj]] = nkoppcells + i
        end
    end
    parent.children .= nkoppcells .+ collect(1:4)
    return nothing
end

function refine!(kopp_grid::KoppGrid, dh::DofHandler, cell::Int)
    nkoppcells = length(kopp_grid.kopp_cells)
    parent = kopp_grid.kopp_cells[cell]
    refinement_level = parent.refinement_level + 1
    for i in 1:4
        parent_secquence = Int[]
        children = Int[0,0,0,0]
        rnt = refined_neighborhood_table[i]
        neighbors = [rnt[j] == 0 ? Ferrite.get_refined_neighbor(kopp_grid,parent.neighbors[j],FaceIndex(i,j)) : nkoppcells + rnt[j] for j in 1:4]
        kc = KoppCell(parent_secquence, refinement_level, cell, parent.root,neighbors, children)
        push!(kopp_grid.kopp_cells, kc)

    end
    for i in 1:4
        neighbors = kopp_grid.kopp_cells[nkoppcells + i].neighbors
        for (jj, j) in enumerate(neighbors)
            j == -1  && continue
            neighbor = kopp_grid.kopp_cells[j]
            neighbor.refinement_level != refinement_level && continue
            neighbor.neighbors[face_face_neighbor_local[jj]] = nkoppcells + i
        end
    end
    parent.children .= nkoppcells .+ collect(1:4)
    # From here starts the dof handling shitshow UwU

    ndofs = Ferrite.ndofs(dh)
    sdh_i = dh.cell_to_subdofhandler[cell] 
    sdh = dh.subdofhandlers[sdh_i]
    ndofs_start = 1
    nn = Int[]
    for interpolation in sdh.field_interpolations
        new_dofs = get_shared_dofs(interpolation)
        nn = copy(new_dofs)
        nbf = getnbasefunctions(interpolation)
        parent_dofs = celldofs(dh, cell)
        nn .= [(x <= nbf) ? parent_dofs[x] : nbf - x - 1 for x in new_dofs]
        for child in 1:4
            neighbors = kopp_grid.kopp_cells[nkoppcells + child].neighbors
            dofs = @view nn[(child-1)*nbf+1:child*nbf]
            for face in boundary_faces[child]
                face_dofs_local = facedof_indices(interpolation)[face]
                neighbors[face] == -1 && continue
                kopp_grid.kopp_cells[neighbors[face]].refinement_level != refinement_level && continue
                o1 = OrientationInfo(faces(kopp_grid.base_grid.cells[parent.root])[face]).flipped
                o2 = OrientationInfo(faces(kopp_grid.base_grid.cells[kopp_grid.kopp_cells[neighbors[face]].root])[face_face_neighbor_local[face]]).flipped
                
                for i in eachindex(face_dofs_local)
                    neighbor_dofs = facedof_indices(interpolation)[face_face_neighbor_local[face]]
                    neighbor_dof = neighbor_dofs[o1 == o2 ? i : length(neighbor_dofs) - i + 1]
                    replace!(nn, dofs[face_dofs_local[i]] => celldofs(dh, neighbors[face])[neighbor_dof])
                end
            end
            for i in eachindex(nn)
                if nn[i] < 0
                    replace!(nn, nn[i] => ndofs + 1)
                    ndofs += 1
                    dh.ndofs[] += 1
                end
            end
        end
        ndofs_start += nbf
    end
    append!(dh.cell_dofs_offset, [length(dh.cell_dofs) + ndofs_per_cell(sdh)*(i-1) + 1 for i in 1:4])
    append!(dh.cell_to_subdofhandler, [dh.cell_to_subdofhandler[cell] for i in 1:4])
    append!(dh.cell_dofs, nn)
    
    return nothing
end