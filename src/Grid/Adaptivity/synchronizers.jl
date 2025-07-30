abstract type AbstractAMRSynchronizer end

struct NullAMRSynchronizer <: AbstractAMRSynchronizer end

struct LTSAMRSynchronizer{V,VCT<:ValuesCache,DHT<:DofHandler} <: AbstractAMRSynchronizer
    values_cache::VCT
    dh::DHT
    celldofs_prev::Vector{Int}
    cell_dofs_offset_prev::Vector{Int}
    cell_to_subdofhandler_prev::Vector{Int}
    cell_matrices::ElasticArray{Float64,3,2,V}
    interface_matrices::ElasticArray{Float64,3,2,V}
    u::Vector{Float64}
    interface_matrix_index::Vector{Int}
end

function sync_amr_refinement_forward!(grid::KoppGrid, sync::LTSAMRSynchronizer, refinement_cache::KoppRefinementCache, n_refined_cells::Int, n_neighborhoods::Int)
    Dim = Ferrite.getspatialdim(grid.base_grid)
    new_length = last(size(sync.cell_matrices))+(2^Dim)*n_refined_cells
    resize!(sync.cell_matrices, (size(sync.cell_matrices)[1], size(sync.cell_matrices)[2], new_length))
    resize!(sync.interface_matrix_index, n_neighborhoods)
    resize!(sync.interface_matrices, (size(sync.interface_matrices)[1], size(sync.interface_matrices)[2], n_neighborhoods ÷ 2))
    # resize!(sync.ansatz_isactive, new_length * sync.dh.subdofhandlers[1].ndofs_per_cell)
    resize!(sync.dh.cell_dofs, new_length * sync.dh.subdofhandlers[1].ndofs_per_cell)
    resize!(sync.dh.cell_dofs_offset, new_length)
    resize!(sync.dh.cell_to_subdofhandler, new_length)

    update_dofs!(grid, refinement_cache, sync.dh, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev)

    sync.interface_matrix_index .= 0
end

function sync_amr_refinement_backward!(sync::LTSAMRSynchronizer)
    copy!(sync.celldofs_prev, sync.dh.cell_dofs)
    copy!(sync.cell_dofs_offset_prev, sync.dh.cell_dofs_offset)
    copy!(sync.cell_to_subdofhandler_prev, sync.dh.cell_to_subdofhandler)

end
