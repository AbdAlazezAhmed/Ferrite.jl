using Ferrite, StaticArrays, ElasticArrays, Tensors

struct KoppCell{Dim, T <: Integer} <: Ferrite.AbstractCell{Ferrite.RefHypercube{Dim}}
    parent::T
    sequence::T
    isleaf::Bool
end

struct ValuesCache{CellValuesT <: CellValues, FacetValuesT <: FacetValues, InterfaceValuesT <: InterfaceValues}
    cell_values::CellValuesT
    facet_values::FacetValuesT
    interface_values::InterfaceValuesT
end

struct KoppGrid{Dim, G <: Ferrite.AbstractGrid{Dim}} <: Ferrite.AbstractGrid{Dim}
    base_grid::G
    kopp_cells::Vector{KoppCell{Dim, Int}}
    kopp_cells_prev::Vector{KoppCell{Dim, Int}}
end

struct KoppRefinementCache{IntT <:Int}
    old_cell_to_new_cell_map::Vector{IntT}
    new_cell_to_old_cell_map::Vector{IntT}
    interfaces_updated_indices::Vector{IntT}
    interfaces_data_updated_indices::Vector{IntT}
    marked_for_refinement::Vector{Bool} #BitVector
    marked_for_coarsening::Vector{Bool} #BitVector
    ncoarseninglevels::IntT
end

include("topology.jl")
include("values.jl")
include("iterators.jl")
include("synchronizers.jl")
include("utils.jl")
include("refinement.jl")
include("coarsening.jl")


function get_refinement_level(cell::KoppCell{Dim, T}) where {Dim, T}
    level::T = 0
    nbits::T = Dim + 1
    mask = (1 << nbits) - 1 # Does this need to be T too?
    maximum_level::T = sizeof(T)*8 ÷ nbits # Maybe use ÷ ?
    for level::T in 0:maximum_level # There should be an easier way
        cell.sequence & (mask << (nbits * level)) == 0 && return level
    end
    return maximum_level
end

function KoppRefinementCache(grid::KoppGrid, topology::KoppTopology)
    n_interfaces = 0
    for _ in InterfaceIterator(grid, topology)
        n_interfaces += 1
    end
    return KoppRefinementCache(
        collect(1:length(topology.root_idx)),
        zeros(Int64,length(topology.root_idx)),
        collect(1:n_interfaces),
        collect(1:n_interfaces),
        zeros(Bool, length(topology.root_idx)),
        zeros(Bool, length(topology.root_idx)),
        # falses(length(topology.root_idx)),
        1)
end



@inline Ferrite.getncells(grid::KoppGrid) = length(grid.kopp_cells)
@inline Ferrite.getcells(grid::KoppGrid, v::Union{Int, Vector{Int}}) = grid.kopp_cells[v]
# @inline getcells(grid::KoppGrid, setname::String) = grid.cells[collect(getcellset(grid,setname))]
@inline Ferrite.getcelltype(grid::KoppGrid) = eltype(grid.kopp_cells)
@inline Ferrite.getcelltype(grid::KoppGrid, i::Int) = typeof(grid.kopp_cells[i])

@inline Ferrite.getrefshape(::KoppCell{Dim}) where Dim = Ferrite.RefHypercube{Dim}

@inline Ferrite.getnnodes(grid::KoppGrid) = length(grid.base_grid.nodes)

function Ferrite.nnodes_per_cell(grid::KoppGrid)
    if !isconcretetype(getcelltype(grid))
        error("There are different celltypes in the `grid`. Use `nnodes_per_cell(grid, cellid::Int)` instead")
    end
    return Ferrite.nnodes(first(grid.base_grid.cells))
end
# TODO: URGENT FIX ME... HOW?
@inline Ferrite.nnodes_per_cell(grid::KoppGrid, i::Int) = Ferrite.nnodes(grid.base_grid.cells[1])

function Ferrite.get_grid(dh::Ferrite.AbstractDofHandler)
    return dh.grid
end

function Ferrite.generate_grid(C::Type{KoppCell{2, T}}, nel::NTuple{2,Int}, left::Vec{2, T2} = Vec((-1.0, -1.0)), right::Vec{2, T2} = Vec((1.0, 1.0))) where {T, T2}
    grid =  generate_grid(Quadrilateral, nel, left, right)
    return KoppGrid(grid)
end

function Ferrite.generate_grid(C::Type{KoppCell{3, T}}, nel::NTuple{3,Int}, left::Vec{3, T2} = Vec((-1.0, -1.0, -1.0)), right::Vec{3, T2} = Vec((1.0, 1.0, 1.0))) where {T, T2}
    grid =  generate_grid(Hexahedron, nel, left, right)
    return KoppGrid(grid)
end

function Ferrite.getneighborhood(topology::KoppTopology, grid::KoppGrid, facetidx::FacetIndex)
    offset = topology.cell_facet_neighbors_offset[facetidx[2], facetidx[1]]
    length = topology.cell_facet_neighbors_length[facetidx[2], facetidx[1]]
    return offset == 0 ? (@view topology.neighbors[0 : -1]) : @view topology.neighbors[offset : offset + length - 1]
end
function Ferrite.getneighborhood(topology::KoppTopology, facetidx::FacetIndex)
    offset = topology.cell_facet_neighbors_offset[facetidx[2], facetidx[1]]
    length = topology.cell_facet_neighbors_length[facetidx[2], facetidx[1]]
    return offset == 0 ? (@view topology.neighbors[0 : -1]) : @view topology.neighbors[offset : offset + length - 1]
end

function KoppGrid(grid::G) where {Dim, G <: Ferrite.AbstractGrid{Dim}}
    cells = KoppCell{Dim, Int}[KoppCell{Dim, Int}(-i, 0, true) for i in 1:getncells(grid)]
    return KoppGrid(grid, cells, copy(cells))
end

function refine!(
        grid::KoppGrid{Dim},
        topology::KoppTopology,
        refinement_cache::KoppRefinementCache,
        sync::AbstractAMRSynchronizer,
        cellset::OrderedSet{CellIndex}
        ) where Dim
    interfaces_dict_prev = _calc_interfaces_dict_prev(grid, topology)

    n_refined_cells = length(cellset)
    new_length = length(grid.kopp_cells)+(2^Dim)*n_refined_cells
    # we need to resize this one early
    _resize_marked_for_refinement!(grid, refinement_cache, cellset)
    refinement_cache.old_cell_to_new_cell_map .= 1:length(refinement_cache.old_cell_to_new_cell_map)
    # Counting how many refined cells and interfaces and calculating the new index
    _update_refinement_cache_isactive!(grid, topology, refinement_cache, cellset)

    _resize_topology!(topology, new_length, Val(Dim))

    zero_topology!(topology)

    n_neighborhoods = count_neighbors_update_indexing!(grid, topology, refinement_cache)
    # sync_amr_refinement_forward!()
    sync_amr_refinement_forward!(grid, sync, refinement_cache, n_refined_cells, n_neighborhoods)

    # Resizing the vectors

    _resize_bunch_of_stuff!(grid, topology, n_neighborhoods, new_length)

    update_cells!(grid, refinement_cache)

    update_root_idx!(grid, topology, refinement_cache)


    # kopp_cache.ansatz_isactive[ndofs_old : end] .= true
    # kopp_cache.interface_matrix_index .= 0


    # # deactivate all parents' shape functions
    update_neighbors!(grid, topology, refinement_cache)
    _calc_interfaces_map(grid, topology, refinement_cache, interfaces_dict_prev)

    # # Refine KoppCache
    # resize!(sync.u, new_length * dh.subdofhandlers[1].ndofs_per_cell)

    # update_koppcache!(grid, refinement_cache, topology, temp_topology, sync, kopp_values, dh, NFacets)

    sync_amr_refinement_backward!(sync, refinement_cache, grid, topology)

    # # # resize refinement cache
    resize!(refinement_cache.old_cell_to_new_cell_map, new_length)
    resize!(refinement_cache.marked_for_refinement, new_length)
    resize!(refinement_cache.marked_for_coarsening, new_length)

    resize!(grid.kopp_cells_prev, length(grid.kopp_cells))
    copy!(grid.kopp_cells_prev, grid.kopp_cells)
    resize!(topology.cell_facet_neighbors_length_prev, size(topology.cell_facet_neighbors_length))
    copy!(topology.cell_facet_neighbors_length_prev, topology.cell_facet_neighbors_length)
    resize!(topology.cell_facet_neighbors_offset_prev, size(topology.cell_facet_neighbors_offset))
    copy!(topology.cell_facet_neighbors_offset_prev, topology.cell_facet_neighbors_offset)
    resize!(topology.neighbors_prev, length(topology.neighbors))
    copy!(topology.neighbors_prev, topology.neighbors)
    resize!(topology.root_idx_prev, length(topology.root_idx))
    copy!(topology.root_idx_prev, topology.root_idx)

    resize!(refinement_cache.interfaces_data_updated_indices, n_neighborhoods ÷ 2)
    refinement_cache.interfaces_data_updated_indices .= 0
    refinement_cache.marked_for_refinement .= false
    refinement_cache.marked_for_coarsening .= false
    return nothing
end

function coarsen!(
        grid::KoppGrid{Dim},
        topology::KoppTopology,
        refinement_cache::KoppRefinementCache,
        sync::AbstractAMRSynchronizer,
        cellset::OrderedSet{CellIndex}
        ) where Dim
    interfaces_dict_prev = _calc_interfaces_dict_prev(grid, topology)
    n_coarsened_cells = length(cellset)
    new_length = length(grid.kopp_cells)-(2^Dim)*n_coarsened_cells
    # we need to resize this one early
    __resize_marked_for_refinement!(grid, refinement_cache, cellset)
    refinement_cache.old_cell_to_new_cell_map .= 1:length(refinement_cache.old_cell_to_new_cell_map)
    # Counting how many refined cells and interfaces and calculating the new index
    __update_refinement_cache_isactive!(grid, topology, refinement_cache, cellset)
    _resize_topology!(topology, new_length, Val(Dim))

    zero_topology!(topology)
    # _resize_bunch_of_stuff!(grid, topology, n_neighborhoods, new_length)

    n_neighborhoods = _count_neighbors_update_indexing!(grid, topology, refinement_cache)
    sync_amr_coarsening_forward!(grid, sync, refinement_cache, n_coarsened_cells, n_neighborhoods)

    resize!(grid.kopp_cells, new_length)
    resize!(topology.root_idx, new_length)
    update_coarsened_cells!(grid, refinement_cache)
    update_coarsened_root_idx!(grid, topology, refinement_cache)




    sync_amr_coarsening_backward!(sync, refinement_cache, grid, topology)

    resize!(refinement_cache.old_cell_to_new_cell_map, new_length)
    resize!(refinement_cache.marked_for_refinement, new_length)
    resize!(refinement_cache.marked_for_coarsening, new_length)

    # update_coarsened_dofs!(grid, refinement_cache, dh, temp_celldofs, temp_cell_dofs_offset, temp_cell_to_subdofhandler)
    # resize!(dh.cell_dofs, new_length * dh.subdofhandlers[1].ndofs_per_cell)
    # resize!(dh.cell_dofs_offset, new_length)
    # resize!(dh.cell_to_subdofhandler, new_length)
    # kopp_cache.interface_matrix_index .= 0

    # # Refine KoppCache
    # resize!(kopp_cache.u, new_length * dh.subdofhandlers[1].ndofs_per_cell)

    # update_coarsened_koppcache!(grid, refinement_cache, topology, temp_topology, kopp_cache, kopp_values, dh, NFacets, Dim)
    # resize!(kopp_cache.interface_matrix_index, n_neighborhoods)
    # resize!(kopp_cache.interface_matrices,size(kopp_cache.interface_matrices)[1], size(kopp_cache.interface_matrices)[2], n_neighborhoods ÷ 2)
    # resize!(kopp_cache.cell_matrices, size(kopp_cache.cell_matrices)[1], size(kopp_cache.cell_matrices)[2], new_length)
    # resize!(kopp_cache.ansatz_isactive, new_length * dh.subdofhandlers[1].ndofs_per_cell)

    resize!(grid.kopp_cells_prev, length(grid.kopp_cells))
    copy!(grid.kopp_cells_prev, grid.kopp_cells)
    resize!(topology.cell_facet_neighbors_length_prev, size(topology.cell_facet_neighbors_length))
    copy!(topology.cell_facet_neighbors_length_prev, topology.cell_facet_neighbors_length)
    resize!(topology.cell_facet_neighbors_offset_prev, size(topology.cell_facet_neighbors_offset))
    copy!(topology.cell_facet_neighbors_offset_prev, topology.cell_facet_neighbors_offset)
    resize!(topology.neighbors_prev, length(topology.neighbors))
    copy!(topology.neighbors_prev, topology.neighbors)
    resize!(topology.root_idx_prev, length(topology.root_idx))
    copy!(topology.root_idx_prev, topology.root_idx)


    resize!(refinement_cache.interfaces_data_updated_indices, n_neighborhoods ÷ 2)
    refinement_cache.interfaces_data_updated_indices .= 0
    # # # # resize refinement cache
    # resize!(refinement_cache.interfaces_updated_indices, length(kopp_cache.interface_matrices) + n_refined_interfaces * (2^(Dim-1) - 1))
    refinement_cache.marked_for_coarsening .= false
    return nothing
end
function Ferrite.getcoordinates!(coords, grid::KoppGrid{Dim}, i) where Dim
    coords .= reinterpret(Vec{Dim, Float64},SVector(get_refined_coords(grid, i)))
    return nothing
end

function to_ferrite_grid(grid::KoppGrid{2})
    cells = Quadrilateral[]
    nodes = Node{2, Float64}[]
    for cc in CellIterator(grid, OrderedSet(1:length(grid.kopp_cells)))
        push!(nodes, Node.(getcoordinates(cc))...)
    end
    unique!(nodes)
    for cc in CellIterator(grid, OrderedSet(1:length(grid.kopp_cells)))
        idx1 = findfirst(x -> x == Node(getcoordinates(cc)[1]), nodes)
        idx2 = findfirst(x -> x == Node(getcoordinates(cc)[2]), nodes)
        idx3 = findfirst(x -> x == Node(getcoordinates(cc)[3]), nodes)
        idx4 = findfirst(x -> x == Node(getcoordinates(cc)[4]), nodes)
        push!(cells, Quadrilateral((idx1, idx2, idx3, idx4)))
    end

    return Ferrite.Grid(cells, nodes)
end
