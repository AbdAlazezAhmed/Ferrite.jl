struct KoppTopology{T <: Integer, V} <: Ferrite.AbstractTopology
    root_idx::Vector{T}
    root_idx_prev::Vector{T}
    cell_facet_neighbors_offset::ElasticMatrix{Int, V}
    cell_facet_neighbors_offset_prev::ElasticMatrix{Int, V}
    cell_facet_neighbors_length::ElasticMatrix{Int, V}
    cell_facet_neighbors_length_prev::ElasticMatrix{Int, V}
    neighbors::Vector{FacetIndex}
    neighbors_prev::Vector{FacetIndex}
    root_facet_orientation_info::Vector{Ferrite.OrientationInfo} # indexed as root_idx + facet_idx - 1
end

function KoppTopology(grid::KoppGrid{Dim}) where Dim
    base_topology = ExclusiveTopology(grid.base_grid)
    root_idx = collect(1:getncells(grid.base_grid))
    neighbors = collect(reinterpret(FacetIndex, reduce(vcat, Ferrite._get_facet_facet_neighborhood(base_topology, Val(Dim)))))
    cell_facet_neighbors_offset = ElasticMatrix(transpose([isempty(neighbor) ? 0 : findfirst(==(FacetIndex(neighbor[1][1], neighbor[1][2])), neighbors) for neighbor in Ferrite._get_facet_facet_neighborhood(base_topology, Val(Dim))]))
    cell_facet_neighbors_length = ElasticMatrix(transpose([isempty(neighbor) ? 0 : 1 for neighbor in Ferrite._get_facet_facet_neighborhood(base_topology, Val(Dim))]))
    root_facet_orientation_info = reduce(vcat, [Ferrite.OrientationInfo(Ferrite.facets(cell)[facet]) for facet in 1:2Dim, cell in getcells(grid.base_grid)])
    return KoppTopology(root_idx, copy(root_idx), cell_facet_neighbors_offset, copy(cell_facet_neighbors_offset), cell_facet_neighbors_length, copy(cell_facet_neighbors_length), neighbors, copy(neighbors), root_facet_orientation_info)
end

@inline function _iterate(topology::KoppTopology, state::NTuple{2, Int} = (1,1))
    state[1] > size(topology.cell_facet_neighbors_offset, 2) && return nothing
    neighborhood = getneighborhood(topology, FacetIndex(state[1], state[2]))
    if state[2] < size(topology.cell_facet_neighbors_offset, 1)
        ret = (FacetIndex(state[1], state[2]), neighborhood)
        state = (state[1], state[2] + 1)
    else
        ret = (FacetIndex(state[1], state[2]), neighborhood)
        state = (state[1] + 1, 1)
    end
    return (ret, state)
end


function Base.iterate(topology::KoppTopology, state::NTuple{2, Int} = (1,1))
    _iterate(topology, state)
end
