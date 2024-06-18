abstract type AbstractAdaptiveGrid{dim} <: AbstractGrid{dim} end
abstract type AbstractAdaptiveCell{refshape <: AbstractRefShape} <: AbstractCell{refshape} end
const refinable_interpolations = [
    Lagrange
    DiscontinuousLagrange
]
function transform_refined(::Lagrange{RefQuadrilateral}, coord::Vec, i::Int)
    displacements = (Vec((-0.5,-0.5)),Vec((0.5,-0.5)),Vec((0.5,0.5)),Vec((-0.5,0.5)))
    return coord/2 + displacements[i]
end

function transform_refined(coord::Vector{<:Vec}, i::Int)
    return [(coord[i]+ coord[j])/2 for j in 1:4]
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
        # @generated function get_boundary_dofs(::$IP)
        #     ref = reference_coordinates($IP())
        #     j = length(ref) + 1
        #     res = Tuple{Int, Int}[]
        #     for i in 1:4
        #         for x in transform_refined.(Ref($IP()), reference_coordinates($IP()), i)
        #             temp = findfirst(u -> x ≈ u, unique(ref))
        #             if temp === nothing
        #                 push!(res,  j)
        #                 j+=1
        #             else
        #                 push!(res,  temp)
        #             end
        #             push!(ref, x)
        #         end
        #     end
        #     return res
        # end
    end
end

mutable struct KoppCell <: AbstractCell{AbstractRefShape{2}}
    secquence::Int # TODO: use two integers instead? alloc == bad
    parent::Int
    children::NTuple{4, Int}
    root::Int
    neighbors::NTuple{4, Pair{Int, Int}}
    isleaf::Bool
    refine::Bool
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

function materialize(kopp_grid::KoppGrid)
    #This is bad, just for viz now
    base_grid = kopp_grid.base_grid
    nodes = copy(base_grid.nodes)
    cells = copy(base_grid.cells)
    for (i, cell) in enumerate(base_grid.cells)
        j = 1
        for leaf in kopp_grid.kopp_roots[i].children
            leaf.isleaf && leaf.parent != -1 || continue
            coords = nodes[[cell.nodes...]]
            parent = leaf
            seq = Int[]
            while true # Not HALAL
                parent_idx = parent.parent
                parent_idx == -1 && break
                push!(seq, parent.secquence)
                parent = kopp_grid.kopp_roots[i].children[parent_idx + 1]
            end
            for k in length(seq):-1:1
                coords = Node.(transform_refined([coord.x for coord in coords], seq[k]))
            end
            l = length(nodes)
            append!(nodes, coords)
            push!(cells, Quadrilateral(ntuple(x -> l + x, 4)))
            j+=1
        end
    end
    Grid(cells, nodes)
    #TODO: check for duplicate nodes
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
            ntuple(j -> 0, 4),
            cell,
            neighborhood,
            true,
            false
        )
        ],ScalarWrapper(0)))
    end
    return Ferrite.KoppGrid(grid, cells)
end
function get_refined_neighbor(grid::G, neighboring_cell::Pair{Int, Int}, current_face::FaceIndex, n::Pair{Int, Int}) where G <: KoppGrid
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

function get_neighbors(grid::KoppGrid, cell::Pair{Int, Int})
    return grid.kopp_roots[cell[1]].children[cell[2] + 1].neighbors
end


function refine!(kopp_grid::KoppGrid, cell_idx::Pair{Int, Int})
    kopp_root = kopp_grid.kopp_roots[cell_idx[1]]
    idx = kopp_root.refinement_index[]
    cell = kopp_root.children[cell_idx[2] + 1]
    @assert cell.isleaf "Cell must be a leaf to be refined"
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
            cell_idx[1],
            ntuple( j -> refined_neighborhood_table[j, i] == 0 ? get_refined_neighbor(kopp_grid, cell.neighbors[j], FaceIndex(i,j), cell_idx[1] => idx + i) : (cell_idx[1] => (idx + refined_neighborhood_table[j,i])) , 4),
            true,
            false
        ))
    end
    cell.children = ntuple(i-> idx + i, 4)
    cell.isleaf = false
    kopp_root.refinement_index[] += 4
end



struct KoppDofHandler{dim,G<:AbstractGrid{dim}} <: AbstractDofHandler
    subdofhandlers::Vector{SubDofHandler{DofHandler{dim, G}}}
    field_names::Vector{Symbol}
    # Dofs for cell i are stored in cell_dofs at the range:
    #     cell_dofs_offset[i]:(cell_dofs_offset[i]+ndofs_per_cell(dh, i)-1)
    cell_dofs::Vector{Vector{Int}}
    cell_dofs_offset::Vector{Vector{Int}}
    cell_to_subdofhandler::Vector{Vector{Int}} # maps cell id -> SubDofHandler id
    closed::ScalarWrapper{Bool}
    grid::G
    ndofs::ScalarWrapper{Int}
end
# function KoppDofHandler(grid::KoppGrid)
#     dh = DofHandler(grid.base_grid)
#     return KoppDofHandler(
#         dh.subdofhandlers,
#         dh.field_names,
#         [Int[] for i in eachindex(grid.kopp_roots)],
#         [[dh.cell_dofs_offset[i]] for i in eachindex(grid.kopp_roots)],
#         [[dh.cell_to_subdofhandler[i]] for i in eachindex(grid.kopp_roots)],
#         dh.closed,
#         grid.base_grid,
#         dh.ndofs
#     )
# end
function KoppDofHandler(grid::KoppGrid, dh::DH) where DH<:DofHandler
    @assert isclosed(dh) "DofHandler must be closed to be converted for adaptivity"
    return KoppDofHandler(
        dh.subdofhandlers,
        dh.field_names,
        # [celldofs(dh, i) for i in eachindex(grid.kopp_roots)],
        [[j for j in eachindex(celldofs(dh, i))] for i in eachindex(grid.kopp_roots)],
        [[1] for i in eachindex(grid.kopp_roots)],
        [[dh.cell_to_subdofhandler[i]] for i in eachindex(grid.kopp_roots)],
        dh.closed,
        dh.grid,
        dh.ndofs
    )
end

function celldofs(dh::KoppDofHandler, cell_idx::Pair{Int, Int})
    sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx[1]][1]]
    start_dof = dh.cell_dofs_offset[cell_idx[1]][cell_idx[2] + 1]
    dh.cell_dofs[cell_idx[1]][start_dof:start_dof + sdh.ndofs_per_cell[]-1]
end

function refine!(kopp_grid::KoppGrid, dh::KoppDofHandler, cell_idx::Pair{Int, Int})
    # Refine the cell and topology
    kopp_root = kopp_grid.kopp_roots[cell_idx[1]]
    idx = kopp_root.refinement_index[]
    cell = kopp_root.children[cell_idx[2] + 1]
    @assert cell.isleaf "Cell must be a leaf to be refined"
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
            cell_idx[1],
            ntuple( j -> refined_neighborhood_table[j, i] == 0 ? get_refined_neighbor(kopp_grid, cell.neighbors[j], FaceIndex(i,j), cell_idx[1] => idx + i) : (cell_idx[1] => (idx + refined_neighborhood_table[j,i])) , 4),
            true,
            false
        ))
        j = 1
        # Refine dofs
        maxdof = maximum(dh.cell_dofs[cell_idx[1]])
        push!(dh.cell_dofs_offset[cell_idx[1]], length(dh.cell_dofs[cell_idx[1]]) + 1)
        sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx[1]][1]]
        for (ip_idx, ip) in enumerate(sdh.field_interpolations)
            n_dofs = getnbasefunctions(ip) * n_components(sdh, ip_idx)
            shared_dofs = get_shared_dofs(ip)[(i-1)*n_dofs + 1 : i*n_dofs]
            for dof in shared_dofs
                push!(dh.cell_dofs[cell_idx[1]], dof > n_dofs ? maxdof + j : celldofs(dh, cell_idx)[dof])
                dof > n_dofs && (j += 1)
            end
        end
    end
    cell.children = ntuple(i-> idx + i, 4)
    cell.isleaf = false
    kopp_root.refinement_index[] += 4

end
