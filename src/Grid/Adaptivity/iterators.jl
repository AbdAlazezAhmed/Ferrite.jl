function Ferrite.CellCache(grid::KoppGrid{dim}, flags::UpdateFlags=Ferrite.UpdateFlags(false, true, false)) where {dim}
    N = Ferrite.nnodes_per_cell(grid, 1) # nodes and coords will be resized in `reinit!`
    nodes = Int[]
    coords = zeros(Vec{dim,Float64}, N)
    return Ferrite.CellCache(flags, grid, -1, nodes, coords, nothing, Int[])
end

function Ferrite.CellCache(dh::DofHandler{dim,<:KoppGrid}, flags::UpdateFlags=Ferrite.UpdateFlags(false, true, true)) where {dim}
    n = Ferrite.ndofs_per_cell(dh.subdofhandlers[1]) # dofs and coords will be resized in `reinit!`
    N = Ferrite.nnodes_per_cell(Ferrite.get_grid(dh), 1)
    nodes = Int[]
    coords = zeros(Vec{dim,Float64}, N)
    celldofs = zeros(Int, n)
    return Ferrite.CellCache(flags, Ferrite.get_grid(dh), -1, nodes, coords, dh, celldofs)
end

function Ferrite.CellCache(sdh::SubDofHandler{<:DofHandler{dim,<:KoppGrid}}, flags::Ferrite.UpdateFlags=UpdateFlags()) where dim
    Tv = Float64
    Ferrite.CellCache(flags, sdh.dh.grid, -1, Int[], Tv[], sdh, Int[])
end


function Ferrite.CellIterator(gridordh::Union{KoppGrid,<:DofHandler{Dim,<:KoppGrid}},
    set::Union{Ferrite.IntegerCollection,Nothing}=nothing,
    flags::UpdateFlags=UpdateFlags(false, true, true)) where Dim
    if set === nothing
        grid = gridordh isa DofHandler ? Ferrite.get_grid(gridordh) : gridordh
        set = 1:getncells(grid)
    end
    if gridordh isa DofHandler
        # TODO: Since the CellCache is resizeable this is not really necessary to check
        #       here, but might be useful to catch slow code paths?
        # _check_same_celltype(get_grid(gridordh), set)
    end
    return Ferrite.CellIterator(CellCache(gridordh, flags), set)
end
function Ferrite.CellIterator(gridordh::Union{KoppGrid,<:DofHandler{Dim,<:KoppGrid}}, flags::UpdateFlags) where Dim
    return Ferrite.CellIterator(gridordh, nothing, flags)
end

function Ferrite.InterfaceIterator(gridordh::Union{Ferrite.AbstractGrid,Ferrite.AbstractDofHandler},
    topology::Ferrite.AbstractTopology=KoppTopology(gridordh isa Ferrite.AbstractGrid ? gridordh : Ferrite.get_grid(gridordh)))
    grid = gridordh isa Ferrite.AbstractGrid ? gridordh : Ferrite.get_grid(gridordh)
    return Ferrite.InterfaceIterator(Ferrite.InterfaceCache(gridordh), grid, topology)
end

function Base.iterate(ii::Ferrite.InterfaceIterator{IC,<:KoppGrid{sdim}}, state::NTuple{3,Int}=(1, 1, 1)) where {sdim,IC}
    ninterfaces = state[3]
    facet_idx = state[1]
    nneighbors = state[2]
    while true
        facet_idx > length(ii.topology.cell_facet_neighbors_offset) && return nothing
        neighborhood_offset = ii.topology.cell_facet_neighbors_offset[(facet_idx-1)%(2*sdim)+1, (facet_idx-1)÷(2*sdim)+1]
        neighborhood_length = ii.topology.cell_facet_neighbors_length[(facet_idx-1)%(2*sdim)+1, (facet_idx-1)÷(2*sdim)+1]
        if neighborhood_offset == 0
            facet_idx = facet_idx + 1
            continue
        end
        neighborhood = @view ii.topology.neighbors[neighborhood_offset:neighborhood_offset+neighborhood_length-1]
        if nneighbors < neighborhood_length
            neighbor = neighborhood[nneighbors]
            facet_idx = facet_idx
            nneighbors = nneighbors + 1
        else
            neighbor = neighborhood[nneighbors]
            facet_idx = facet_idx + 1
            nneighbors = 1
        end
        _facet_idx = nneighbors == 1 ? facet_idx - 1 : facet_idx
        _nneighbors = nneighbors == 1 ? nneighbors : nneighbors - 1
        cell_idx = (_facet_idx - 1) ÷ (2 * sdim) + 1
        facet_a = (_facet_idx - 1) % (2 * sdim) + 1
        cell = ii.grid.kopp_cells[cell_idx]
        cell_refinement_level = get_refinement_level(cell)
        neighbor_refinement_level = get_refinement_level(ii.grid.kopp_cells[neighbor[1]])
        if ((cell_refinement_level == neighbor_refinement_level) && (neighbor[1] > cell_idx)) || (cell_refinement_level < neighbor_refinement_level)
            reinit!(ii.cache, FacetIndex(cell_idx, facet_a), neighbor)
            ninterfaces += 1
            return (facet_idx, nneighbors, ninterfaces), (facet_idx, nneighbors, ninterfaces)
        else
            continue
        end
    end
    return nothing
end
