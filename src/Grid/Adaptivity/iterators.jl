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
        set = findall(cell -> cell.isleaf, grid.kopp_cells)
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

@inline function interface_topology_iterate(topology::KoppTopology, state::NTuple{2, Int} = (1,1))
    state[1] > size(topology.cell_facet_neighbors_offset, 2) && return ((FacetIndex(0, 0), FacetIndex(0, 0)), (-1, -1))
    if topology.cell_facet_neighbors_length[state[2], state[1]] != 1
        if state[2] < size(topology.cell_facet_neighbors_offset, 1)
            ret = (FacetIndex(state[1], state[2]), FacetIndex(0,0))
            state = (state[1], state[2] + 1)
        else
            ret = (FacetIndex(state[1], state[2]), FacetIndex(0,0))
            state = (state[1] + 1, 1)
        end
    else
        neighborhood = topology.neighbors[topology.cell_facet_neighbors_offset[state[2], state[1]]]
        if state[2] < size(topology.cell_facet_neighbors_offset, 1)
            ret = (FacetIndex(state[1], state[2]), neighborhood)
            state = (state[1], state[2] + 1)
        else
            ret = (FacetIndex(state[1], state[2]), neighborhood)
            state = (state[1] + 1, 1)
        end
    end

    return (ret, state)
end

function Base.iterate(ii::Ferrite.InterfaceIterator{IC,<:KoppGrid{sdim}}, state = (1,1)) where {sdim,IC}
        while true
            it = interface_topology_iterate(ii.topology, state)
            (facet_a, neighbors), state = it
            state == (-1, -1) && return nothing
            if neighbors == FacetIndex(0,0)
                continue
            end
            facet_b = neighbors
            get_refinement_level(ii.grid.kopp_cells[facet_a[1]]) == get_refinement_level(ii.grid.kopp_cells[facet_b[1]]) && facet_a[1] > facet_b[1] && continue
            reinit!(ii.cache, facet_a, facet_b)
            return (ii.cache, state)
        end

    return
end
