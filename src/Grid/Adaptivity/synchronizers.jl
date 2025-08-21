using Unrolled

abstract type AbstractAMRSynchronizer end

struct NullAMRSynchronizer <: AbstractAMRSynchronizer end

abstract type AMRDataContainer end

abstract type AMRCellData{C} <: AMRDataContainer end

struct CellMassMatrix{C} <: AMRCellData{C}
    data::C
end

struct CellStiffnessMatrix{C} <: AMRCellData{C}
    data::C
end

abstract type AMRInterfaceData{C} <: AMRDataContainer end

struct InterfaceStiffnessMatrix{C} <: AMRInterfaceData{C}
    data::C
end

abstract type AMRdofwiseVector{C} <: AMRDataContainer end

struct SolutionVector{C} <: AMRdofwiseVector{C}
    data::C
end

struct TimeVector{C} <: AMRCellData{C}
    data::C
end

struct ErrorVector{C} <: AMRCellData{C}
    data::C
end

function Base.resize!(data_store::ErrorVector{<:ElasticArray}, data_store_prev::ErrorVector{<:ElasticArray}, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store_prev.data, size(data_store.data))
    data_store_prev.data .= data_store.data
    Base.resize!(data_store.data, ncells)
    data_store.data .= zero(eltype(data_store.data))
end

struct CellNeighboringDofsMap
    offsets::Vector{Int}
    dofs::Vector{Int}
end

function Base.resize!(data_store::CellNeighboringDofsMap, data_store_prev::CellNeighboringDofsMap, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store_prev.offsets, size(data_store.offsets))
    Base.resize!(data_store_prev.dofs, size(data_store.dofs))
    data_store_prev.dofs .= data_store.dofs
    data_store_prev.offsets .= data_store.offsets
    Base.resize!(data_store.dofs, ncells * ndofs_per_cell + ninterfaces * 2 * ndofs_per_cell)
    Base.resize!(data_store.offsets, ncells)
    data_store.offsets .= zero(eltype(data_store.offsets))
    data_store.dofs .= zero(eltype(data_store.dofs))
end

struct AssembledStiffnessMatrix{C <: ElasticArray}
    data::C
end

function update!(data_store::AssembledStiffnessMatrix{<:ElasticArray}, data_store_prev::AssembledStiffnessMatrix{<:ElasticArray}, grid, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int}, values_cache, sync)
    interface_index = 1
    dofs_map = sync.data_stores[5]
    K_assembled = sync.data_stores[4]
    ndofs_cell =
    dofs_temp_storage = zeros(Int, ndofs_cell)
    dofs_temp_storage2 = zeros(Int, 2*ndofs_cell)
    @inbounds @views for ic in InterfaceIterator(grid, topology)
        cell_a = cellid(ic.a)
        cell_b = cellid(ic.b)
        K[cell_a][1:ndofs_cell,1:ndofs_cell] .+= sync.data_stores[3].data[1:ndofs_cell,1:ndofs,interface_index]
        celldofs!(dofs_temp_storage, dh, cell_a)
        dofs_map[cell_a][1:ndofs_cell] .= dofs_temp_storage

        neighbor_idx = 1
        found = false
        @inbounds for facet in 1:4
            if found
                break
            end
            neighborhood = getneighborhood(topology, FacetIndex(cell_a, facet))

            @inbounds for (i, neighbor) in enumerate(neighborhood)
                if neighbor[1] == cell_b
                    dof_start = ndofs_cell + 1 + (neighbor_idx-1)*ndofs_cell
                    dof_end = ndofs_cell + neighbor_idx*ndofs_cell
                    K[cell_a][1:ndofs_cell,dof_start:dof_end] += sync.data_stores[3].data[1:ndofs_cell, ndofs_cell + 1:end, interface_index]
                    celldofs!(dofs_temp_storage, dh, neighbor[1])
                    dofs_map[cell_a][dof_start:dof_end] .= dofs_temp_storage
                    found = true
                    break
                end
                neighbor_idx += 1
            end
        end

        K[cell_b][1:ndofs_cell,1:ndofs_cell] .+= sync.data_stores[3].data[ndofs_cell+1:end, ndofs_cell+1:end,interface_index]
        celldofs!(dofs_temp_storage, dh, cell_b)
        dofs_map[cell_b][1:ndofs_cell] .= dofs_temp_storage

        neighbor_idx = 1
        global found = false
        for facet in 1:4
            if found
                break
            end
            neighborhood = getneighborhood(topology, FacetIndex(cell_b, facet))
            for (i, neighbor) in enumerate(neighborhood)
                if neighbor[1] == cell_a
                    dof_start = ndofs_cell + 1 + (neighbor_idx-1)*ndofs_cell
                    dof_end = ndofs_cell + neighbor_idx*ndofs_cell
                    K[cell_b][1:ndofs_cell,dof_start:dof_end] += sync.data_stores[3].data[ndofs_cell + 1:end, 1:ndofs_cell, interface_index]
                    celldofs!(dofs_temp_storage, dh,  neighbor[1])

                    dofs_map[cell_b][dof_start:dof_end] .= dofs_temp_storage
                    found = true
                    break
                end
                neighbor_idx += 1
            end
        end
        interface_index += 1
        reinit!(interfacevalues, ic, topology)
        celldofs!(( dofs_temp_storage2[1:ndofs_cell]), dh, cell_a)
        celldofs!(( dofs_temp_storage2[ndofs_cell + 1:end]), dh, cell_b)

        estimate_kelly_interface!(Float64, error_vector, u[dofs_temp_storage2], ic, sync.values_cache.interface_values)
    end
end


function Base.resize!(data_store::AssembledStiffnessMatrix, data_store_prev::AssembledStiffnessMatrix, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store_prev.data, size(data_store.data))
    data_store_prev.data .= data_store.data
    Base.resize!(data_store.data, size(data_store.data)[1], ncells * ndofs_per_cell + ninterfaces * 2 * ndofs_per_cell)
    data_store.data .= zero(eltype(data_store.data))
end

function assemble_element_matrix!(K_cell_matrices, kopp_values, cell_idx)
    error("FUCK MICROSOFT")
end

function Base.resize!(data_store::TimeVector, data_store_prev::TimeVector, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store_prev.data, length(data_store.data))
    data_store_prev.data .= data_store.data
    Base.resize!(data_store.data, ncells)
    # data_store.data .= zero(eltype(data_store.data))
end

function Base.resize!(data_store::AMRdofwiseVector, data_store_prev::AMRdofwiseVector, ncells, ninterfaces, ndofs)
    Base.resize!(data_store_prev.data, length(data_store.data))
    data_store_prev.data .= data_store.data
    Base.resize!(data_store.data, ndofs)
    # data_store.data .= zero(eltype(data_store.data))
end

function Base.resize!(data_store::AMRCellData{<:ElasticArray}, data_store_prev::AMRCellData{<:ElasticArray}, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store_prev.data, size(data_store.data))
    data_store_prev.data .= data_store.data
    Base.resize!(data_store.data, size(data_store.data)[1], size(data_store.data)[2], ncells)
    data_store.data .= zero(eltype(data_store.data))
end

function Base.resize!(data_store::AMRInterfaceData{<:ElasticArray}, data_store_prev::AMRInterfaceData{<:ElasticArray}, ncells, ninterfaces, ndofs_per_cell)
    Base.resize!(data_store_prev.data, size(data_store.data))
    data_store_prev.data .= data_store.data
    Base.resize!(data_store.data, size(data_store.data)[1], size(data_store.data)[2], ninterfaces)
    data_store.data .= zero(eltype(data_store.data))
end

function update!(data_store::TimeVector, data_store_prev::TimeVector, grid::KoppGrid{Dim}, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int}, values_cache) where Dim
    @views @inbounds for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        new == 0 && continue
        data_store.data[new] = data_store_prev.data[old]
        if refinement_cache.marked_for_refinement[new]
            data_store.data[new+1:new+2^Dim] .= data_store.data[new]
        end
    end
end

function update!(data_store::AMRdofwiseVector, data_store_prev::AMRdofwiseVector, grid::KoppGrid{Dim}, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int}, values_cache) where Dim
    data_store.data .= 0.0
    sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[1]]
    empty!(sdh.cellset)
    @inbounds  @views for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        new == 0 && continue
        new_offset = dh.cell_dofs_offset[new]
        new_sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[new]]
        # TODO: do only once
        new_dofs = dh.cell_dofs[new_offset:new_offset+new_sdh.ndofs_per_cell-1]
        old_offset = cell_dofs_offset_prev[old]
        old_sdh = dh.subdofhandlers[cell_to_subdofhandler_prev[old]]
        old_dofs = celldofs_prev[old_offset:old_offset+old_sdh.ndofs_per_cell-1]
        data_store.data[new_dofs] .= data_store_prev.data[old_dofs]
        # if refinement_cache.marked_for_coarsening[new]
        #     delete!(sdh.cellset, old)
        #     delete!(sdh.cellset, old+1)
        #     delete!(sdh.cellset, old+2)
        #     delete!(sdh.cellset, old+3)
        #     delete!(sdh.cellset, old+4)
        # end
        push!(sdh.cellset, new)
    end

    @inbounds @views for cell_idx in 1:length(grid.kopp_cells)
        cell = getcells(grid, cell_idx)
        # Calculate new interfaces
        if refinement_cache.marked_for_refinement[cell_idx]
            parent_offset = dh.cell_dofs_offset[cell_idx]
            parent_sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cell_idx]]
            parent_dofs = dh.cell_dofs[parent_offset:parent_offset+parent_sdh.ndofs_per_cell-1]
            for new_seq in 1:2^Dim
                child_idx = cell_idx + new_seq
                push!(sdh.cellset, child_idx)
                # DofHandler add new dofs and offsets
                # For DG only
                offset = dh.cell_dofs_offset[child_idx]
                sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[child_idx]]
                dofs = dh.cell_dofs[offset:offset+sdh.ndofs_per_cell-1]
                # for dof in dofs
                #     @info new_seq dof data_store.data[dof]
                # end
                interpolate_solution!(Ferrite.RefHypercube{Dim}, data_store.data, dofs, parent_dofs, new_seq)
            end
        end
        if refinement_cache.marked_for_coarsening[cell_idx]
        end
    end
end

function update!(data_store::AMRCellData{<:ElasticArray}, data_store_prev::AMRCellData{<:ElasticArray}, grid::KoppGrid{Dim}, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int}, values_cache) where Dim
    @views @inbounds for (old, new) in enumerate(refinement_cache.old_cell_to_new_cell_map)
        new == 0 && continue
        data_store.data[:, :, new] .= data_store_prev.data[:, :, old]
        if refinement_cache.marked_for_refinement[new]
            data_store.data[:, :, new+1:new+2^Dim] .= 0.0
        end
    end
    @inbounds for cell_cache in Ferrite.CellIterator(grid)
        cell_idx = cell_cache.cellid
        old = refinement_cache.new_cell_to_old_cell_map[cell_idx]
        if iszero(old)
            reinit!(values_cache.cell_values, cell_cache)
            assemble_element_matrix!(data_store, values_cache, cell_idx)
        else
            new_data = @view data_store.data[:, :, cell_idx]
            new_data .= @view data_store_prev.data[:, :, old]
        end
    end
end

function update!(data_store::AMRInterfaceData{<:ElasticArray}, data_store_prev::AMRInterfaceData{<:ElasticArray}, grid, topology, dh, refinement_cache, celldofs_prev::Vector{Int}, cell_dofs_offset_prev::Vector{Int}, cell_to_subdofhandler_prev::Vector{Int}, values_cache)
    interface_index = 1
    data_store.data .= 0.0
    @inbounds for (old, new) in enumerate(refinement_cache.interfaces_data_updated_indices)
        iszero(new) && continue
        new_data = @view data_store.data[:, :, new]
        new_data .= @view data_store_prev.data[:, :, old]
    end
    @inbounds for interface_cache in Ferrite.InterfaceIterator(grid, topology)
        data = @view data_store.data[:, :, interface_index]
        if all(iszero, data)
            reinit!(values_cache.interface_values, interface_cache, topology)
            assemble_element_matrix!(data_store, values_cache, interface_index)
        end
        interface_index += 1
    end
end

struct LTSAMRSynchronizer{Order,DataStores<:Tuple,VCT<:ValuesCache,DHT<:DofHandler} <: AbstractAMRSynchronizer
    # Values
    values_cache::VCT
    # DofHandler
    dh::DHT
    celldofs_prev::Vector{Int}
    cell_dofs_offset_prev::Vector{Int}
    cell_to_subdofhandler_prev::Vector{Int}
    # Where matrices and vectors are stored
    data_stores_prev::DataStores # Tuple of matrices and/or vectors wrapped in AMRDataContainer
    data_stores::DataStores # Tuple of matrices and/or vectors wrapped in AMRDataContainer
    # Indexing
    interface_matrix_index::Vector{Int}
    # scalars
    Δt::Float64
end

function sync_amr_refinement_forward!(grid::KoppGrid, sync::LTSAMRSynchronizer{Order}, refinement_cache::KoppRefinementCache, n_refined_cells::Int, n_neighborhoods::Int) where Order
    Dim = Ferrite.getspatialdim(grid.base_grid)
    new_length = Ferrite.getncells(grid) + (2^Dim) * n_refined_cells
    n = length(sync.data_stores)
    Unroll.@unroll for i in 1:5
        Base.resize!(sync.data_stores[i], sync.data_stores_prev[i], new_length, n_neighborhoods ÷ 2, maximum(sync.dh.cell_dofs) + n_refined_cells * (2^Dim - 1) * sync.dh.subdofhandlers[1].ndofs_per_cell)
    end


    resize!(sync.celldofs_prev, length(sync.dh.cell_dofs))
    resize!(sync.cell_dofs_offset_prev, length(sync.dh.cell_dofs_offset))
    resize!(sync.cell_to_subdofhandler_prev, length(sync.dh.cell_to_subdofhandler))

    copy!(sync.celldofs_prev, (sync.dh.cell_dofs))
    copy!(sync.cell_dofs_offset_prev, (sync.dh.cell_dofs_offset))
    copy!(sync.cell_to_subdofhandler_prev, (sync.dh.cell_to_subdofhandler))

    resize!(sync.dh.cell_dofs, new_length * sync.dh.subdofhandlers[1].ndofs_per_cell)
    resize!(sync.dh.cell_dofs_offset, new_length)
    resize!(sync.dh.cell_to_subdofhandler, new_length)

    update_dofs!(grid, refinement_cache, sync.dh, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev, Val(Order))

    sync.interface_matrix_index .= 0
end

function sync_amr_coarsening_forward!(grid::KoppGrid, sync::LTSAMRSynchronizer{Order}, refinement_cache::KoppRefinementCache, n_coarsened_cells::Int, n_neighborhoods::Int) where Order
    Dim = Ferrite.getspatialdim(grid.base_grid)
    new_length = Ferrite.getncells(grid) - (2^Dim) * n_coarsened_cells
    Unroll.@unroll for i in 1:5
        Base.resize!(sync.data_stores[i], sync.data_stores_prev[i], new_length, n_neighborhoods ÷ 2, maximum(sync.dh.cell_dofs) - n_coarsened_cells * (2^Dim - 1) * sync.dh.subdofhandlers[1].ndofs_per_cell)
    end
    resize!(sync.celldofs_prev, length(sync.dh.cell_dofs))
    resize!(sync.cell_dofs_offset_prev, length(sync.dh.cell_dofs_offset))
    resize!(sync.cell_to_subdofhandler_prev, length(sync.dh.cell_to_subdofhandler))

    copy!(sync.celldofs_prev, (sync.dh.cell_dofs))
    copy!(sync.cell_dofs_offset_prev, (sync.dh.cell_dofs_offset))
    copy!(sync.cell_to_subdofhandler_prev, (sync.dh.cell_to_subdofhandler))

    resize!(sync.dh.cell_dofs, new_length * sync.dh.subdofhandlers[1].ndofs_per_cell)
    resize!(sync.dh.cell_dofs_offset, new_length)
    resize!(sync.dh.cell_to_subdofhandler, new_length)

    update_dofs!(grid, refinement_cache, sync.dh, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev, Val(Order))

    sync.interface_matrix_index .= 0
end

function sync_amr_refinement_backward!(sync::LTSAMRSynchronizer, refinement_cache::KoppRefinementCache, grid, topology)
    Unroll.@unroll for i in 1:5
        update!(sync.data_stores[i], sync.data_stores_prev[i], grid, topology, sync.dh, refinement_cache, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev, sync.values_cache)
    end
    copy!(sync.celldofs_prev, sync.dh.cell_dofs)
    copy!(sync.cell_dofs_offset_prev, sync.dh.cell_dofs_offset)
    copy!(sync.cell_to_subdofhandler_prev, sync.dh.cell_to_subdofhandler)

end

function sync_amr_coarsening_backward!(sync::LTSAMRSynchronizer, refinement_cache::KoppRefinementCache, grid, topology)
    Unroll.@unroll for i in 1:5
        update!(sync.data_stores[i], sync.data_stores_prev[i], grid, topology, sync.dh, refinement_cache, sync.celldofs_prev, sync.cell_dofs_offset_prev, sync.cell_to_subdofhandler_prev, sync.values_cache)
    end

    copy!(sync.celldofs_prev, sync.dh.cell_dofs)
    copy!(sync.cell_dofs_offset_prev, sync.dh.cell_dofs_offset)
    copy!(sync.cell_to_subdofhandler_prev, sync.dh.cell_to_subdofhandler)

end

function assemble_element_matrix!(K::CellMassMatrix, kopp_values::ValuesCache, cell_idx)
    Ke = @view K.data[:, :, cell_idx]
    cv = kopp_values.cell_values
    n_basefuncs = getnbasefunctions(cv)
    for q_point in 1:getnquadpoints(cv.qr)
        dΩ = getdetJdV(cv, q_point)
        for i in 1:n_basefuncs
            δu = Ferrite.shape_value(cv, q_point, i)
            for j in 1:n_basefuncs
                u = Ferrite.shape_value(cv, q_point, j)
                Ke[i, j] += (δu ⋅ u) * dΩ
            end
        end
    end
    return nothing
end
function assemble_element_matrix!(K::CellStiffnessMatrix, kopp_values::ValuesCache, cell_idx)
    Ke = @view K.data[:, :, cell_idx]
    cv = kopp_values.cell_values
    n_basefuncs = getnbasefunctions(cv)
    for q_point in 1:getnquadpoints(cv)
        ## Get the quadrature weight
        dΩ = getdetJdV(cv, q_point)
        ## Loop over test shape functions
        for i in 1:n_basefuncs
            ∇δu = shape_gradient(cv, q_point, i)
            ## Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cv, q_point, j)
                ## Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return nothing
end
function assemble_element_matrix!(K::InterfaceStiffnessMatrix, kopp_values::ValuesCache, interface_index, μ::Float64=8.)
    Ki = @view K.data[:, :, interface_index]
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
                Ki[i, j] += -(δu_jump ⋅ ∇u_avg + ∇δu_avg ⋅ u_jump) * dΓ + μ * (δu_jump ⋅ u_jump) * dΓ
            end
        end
    end
end


function LTSAMRSynchronizer(grid::KoppGrid, dh::Ferrite.AbstractDofHandler, kopp_values::ValuesCache, refinement_cache::KoppRefinementCache, topology::KoppTopology{NFacets}, Δt=0.1, _u=zeros(ndofs(dh))) where NFacets
    @assert Ferrite.isclosed(dh) "DofHandler must be closed"
    ndofs_cell = ndofs_per_cell(dh)
    # M
    M_cell_matrices = CellMassMatrix(ElasticArray(zeros(Float64, ndofs_cell, ndofs_cell, length(refinement_cache.old_cell_to_new_cell_map))))
    for cell_cache in Ferrite.CellIterator(grid)
        reinit!(kopp_values.cell_values, cell_cache)
        cell_idx = cell_cache.cellid
        assemble_element_matrix!(M_cell_matrices, kopp_values, cell_idx)
    end
    # K_cell_matrices
    K_cell_matrices = CellStiffnessMatrix(ElasticArray(zeros(Float64, ndofs_cell, ndofs_cell, length(refinement_cache.old_cell_to_new_cell_map))))
    for cell_cache in Ferrite.CellIterator(grid)
        reinit!(kopp_values.cell_values, cell_cache)
        cell_idx = cell_cache.cellid
        assemble_element_matrix!(K_cell_matrices, kopp_values, cell_idx)
    end
    K_interface_matrices = InterfaceStiffnessMatrix(ElasticArray(zeros(Float64, 2 * ndofs_cell, 2 * ndofs_cell, length(topology.neighbors) ÷ 2)))
    # Construct interface matrix indices
    interface_matrix_indices = Vector{Int}(undef, count(topology.cell_facet_neighbors_length .!= 0))
    interface_index = 1
    for interface_cache in Ferrite.InterfaceIterator(grid, topology)
        reinit!(kopp_values.interface_values, interface_cache, topology)
        assemble_element_matrix!(K_interface_matrices, kopp_values, interface_index)
        interface_index += 1
    end

    # N_cell_vectors = zero(u)
    # new_offset = 1
    # for (idx, offset) in pairs(IndexCartesian(), topology.cell_facet_neighbors_offset)
    #     offset == 0 && continue
    #     neighbor_cell = topology.neighbors[offset][1]
    #     neighbor_facet = topology.neighbors[offset][2]
    #     current_cell = idx[2]
    #     if(neighbor_cell < current_cell)
    #         interface_matrix_indices[offset] = interface_matrix_indices[topology.cell_facet_neighbors_offset[neighbor_facet, neighbor_cell]]
    #         continue
    #     end
    #     interface_matrix_indices[offset] = new_offset
    #     new_offset += 1
    # end
    u = SolutionVector(_u)
    t = TimeVector(zeros(getncells(grid)))
    ansatz_isactive = trues(ndofs(dh))
    datastors = (M_cell_matrices, K_cell_matrices, K_interface_matrices, u, t)
    return LTSAMRSynchronizer{
        Ferrite.getorder(dh.subdofhandlers[dh.cell_to_subdofhandler[1]].field_interpolations[1]),
        typeof(datastors),
        typeof(kopp_values),
        typeof(dh)}(kopp_values,
        dh,
        copy(dh.cell_dofs),
        copy(dh.cell_dofs_offset),
        copy(dh.cell_to_subdofhandler),
        datastors,
        (deepcopy(M_cell_matrices), deepcopy(K_cell_matrices), deepcopy(K_interface_matrices), deepcopy(u), deepcopy(t)),
        interface_matrix_indices,
        Δt)
end

struct KoppSynchronizer end
