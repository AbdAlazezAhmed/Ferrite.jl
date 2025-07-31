using Ferrite, StaticArrays, ElasticArrays, Tensors

struct KoppCell{Dim, T <: Integer} <: Ferrite.AbstractCell{Ferrite.RefHypercube{Dim}}
    parent::T
    sequence::T
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
    children_updated_indices::Vector{IntT}
    interfaces_updated_indices::Vector{IntT}
    marked_for_refinement::BitVector
    marked_for_coarsening::BitVector
    ncoarseninglevels::IntT
end

include("topology.jl")
include("values.jl")
include("iterators.jl")
include("synchronizers.jl")
include("utils.jl")


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
    interfaces = zeros(Int, length(topology.neighbors))
    n_interfaces = 0
    for (i, offset) in pairs(IndexCartesian(), topology.cell_facet_neighbors_offset)
            offset == 0 && continue
            neighbor = topology.neighbors[offset]
            # Refine interface when the current cell is lower index than the neighbor only to avoid race conditions
            get_refinement_level(grid.kopp_cells[neighbor.idx[1]]) == get_refinement_level(grid.kopp_cells[i[2]]) && i[2] > neighbor.idx[1] && continue
            n_interfaces += 1
            interfaces[offset] = n_interfaces
            interfaces[topology.cell_facet_neighbors_offset[neighbor.idx[2], neighbor.idx[1]]] = n_interfaces
    end
    return KoppRefinementCache(
        collect(1:length(topology.root_idx)),
        interfaces,
        falses(length(topology.root_idx)),
        falses(length(topology.root_idx)),
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

function Ferrite.generate_grid(C::Type{KoppCell{2, T}}, nel::NTuple{2,Int}) where {T}
    grid =  generate_grid(Quadrilateral, nel)
    return KoppGrid(grid)
end

function Ferrite.generate_grid(C::Type{KoppCell{3, T}}, nel::NTuple{3,Int}) where {T}
    grid =  generate_grid(Hexahedron, nel)
    return KoppGrid(grid)
end


function Ferrite.getneighborhood(topology::KoppTopology, grid::KoppGrid, facetidx::FacetIndex)
    offset = topology.cell_facet_neighbors_offset[facetidx[2], facetidx[1]]
    length = topology.cell_facet_neighbors_length[facetidx[2], facetidx[1]]
    return offset == 0 ? nothing : @view topology.neighbors[offset : offset + length - 1]
end

function KoppGrid(grid::G) where {Dim, G <: Ferrite.AbstractGrid{Dim}}
    cells = KoppCell{Dim, Int}[KoppCell{Dim, Int}(-i, 0) for i in 1:getncells(grid)]
    return KoppGrid(grid, cells, copy(cells))
end

function refine!(
        grid::KoppGrid{Dim},
        topology::KoppTopology,
        refinement_cache::KoppRefinementCache,
        sync::AbstractAMRSynchronizer,
        cellset::Set{CellIndex}
        ) where Dim
    n_refined_cells = length(cellset)
    new_length = length(grid.kopp_cells)+(2^Dim)*n_refined_cells
    # we need to resize this one early
    _resize_marked_for_refinement!(grid, refinement_cache, cellset)
    n_refined_interfaces = 0
    NFacets = 2*Dim
    refinement_cache.children_updated_indices .= 1:length(refinement_cache.children_updated_indices)
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

    # # Refine KoppCache
    # resize!(sync.u, new_length * dh.subdofhandlers[1].ndofs_per_cell)

    # update_koppcache!(grid, refinement_cache, topology, temp_topology, sync, kopp_values, dh, NFacets)

    # # # resize refinement cache
    resize!(refinement_cache.children_updated_indices, new_length)
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
    refinement_cache.marked_for_refinement .= false
    refinement_cache.marked_for_coarsening .= false
    sync_amr_refinement_backward!(sync)
    return nothing
end

function coarsen!(
        grid::KoppGrid{Dim},
        topology::KoppTopology,
        refinement_cache::KoppRefinementCache,
        kopp_cache::AbstractAMRSynchronizer,
        cellset::Set{CellIndex},
        dh::DofHandler) where Dim
    n_coarsened_cells = length(cellset)
    new_length = length(grid.kopp_cells)-(2^Dim)*n_coarsened_cells
    # we need to resize this one early
    __resize_marked_for_refinement!(grid, refinement_cache, cellset)
    # n_refined_interfaces = 0
    NFacets = 2*Dim
    refinement_cache.children_updated_indices .= 1:length(refinement_cache.children_updated_indices)
    # Counting how many refined cells and interfaces and calculating the new index
    __update_refinement_cache_isactive!(grid, topology, refinement_cache, kopp_cache, cellset, dh)

    temp_grid = deepcopy(grid)
    _resize_topology!(topology, new_length, Val(Dim))

    zero_topology!(topology)

    n_neighborhoods = update_topology!(grid, topology, refinement_cache)

    # # Resizing the vectors
    temp_celldofs = copy(dh.cell_dofs)
    temp_cell_dofs_offset = copy(dh.cell_dofs_offset)
    temp_cell_to_subdofhandler = copy(dh.cell_to_subdofhandler)

    ndofs_old = dh.ndofs

    resize!(topology.neighbors, n_neighborhoods)
   resize!(topology.root_idx, new_length)

    update_coarsened_cells!(grid, refinement_cache)
    resize!(grid.kopp_cells, new_length)

    update_coarsened_root_idx!(grid, topology, refinement_cache)

    update_coarsened_dofs!(grid, refinement_cache, dh, temp_celldofs, temp_cell_dofs_offset, temp_cell_to_subdofhandler)
    resize!(dh.cell_dofs, new_length * dh.subdofhandlers[1].ndofs_per_cell)
    resize!(dh.cell_dofs_offset, new_length)
    resize!(dh.cell_to_subdofhandler, new_length)
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

    # # # # resize refinement cache
    resize!(refinement_cache.children_updated_indices, new_length)
    resize!(refinement_cache.marked_for_coarsening, new_length)
    resize!(refinement_cache.marked_for_refinement, new_length)
    # resize!(refinement_cache.interfaces_updated_indices, length(kopp_cache.interface_matrices) + n_refined_interfaces * (2^(Dim-1) - 1))
    refinement_cache.marked_for_coarsening .= false
    return nothing
end
function Ferrite.getcoordinates!(coords, grid::KoppGrid{Dim}, i) where Dim
    coords .= reinterpret(Vec{Dim, Float64},SVector(get_refined_coords(grid, i)))
    return nothing
end

# # # 2 allocs for J?
function assemble_element_matrix!(Ke, kopp_values::ValuesCache)
    cv = kopp_values.cell_values
    n_basefuncs = getnbasefunctions(cv)
    for q_point in 1:getnquadpoints(cv.qr)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs
            δu  = Ferrite.shape_value(cv, q_point, i)
            for j in 1:n_basefuncs
                u = Ferrite.shape_value(cv, q_point, i)
                Ke[i, j] += (δu ⋅ u) * dΩ
            end
        end
    end
    return nothing
end


function assemble_interface_matrix!(Ki, kopp_values::ValuesCache, μ::Float64 = 10.)
    iv = kopp_values.interface_values
    for q_point in 1:getnquadpoints(iv)
        # Get the normal to facet A
        normal = getnormal(iv, q_point)
        # Get the quadrature weight
        dΓ = getdetJdV(iv, q_point)
        # Loop over test shape functions
        for i in 1:getnbasefunctions(iv)
            # Multiply the jump by the negative normal to get the definition from the theory section.
            δu_jump = shape_value_jump(iv, q_point, i) * (-normal)
            ∇δu_avg = shape_gradient_average(iv, q_point, i)
            # Loop over trial shape functions
            for j in 1:getnbasefunctions(iv)
                # Multiply the jump by the negative normal to get the definition from the theory section.
                u_jump = shape_value_jump(iv, q_point, j) * (-normal)
                ∇u_avg = shape_gradient_average(iv, q_point, j)
                # Add contribution to Ki
                Ki[i, j] += -(δu_jump ⋅ ∇u_avg + ∇δu_avg ⋅ u_jump)*dΓ +  μ * (δu_jump ⋅ u_jump) * dΓ
            end
        end
    end
end
