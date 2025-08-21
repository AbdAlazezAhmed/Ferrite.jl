function assemble_element!(Me::Matrix, Ke::Matrix, fe::Vector, cellvalues::CellValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    ## Reset to 0
    fill!(Ke, 0)
    fill!(Me, 0)
    fill!(fe, 0)
    ## Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        ## Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        ## Loop over test shape functions
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            ## Add contribution to fe
            #fe[i] += δu * dΩ
            ## Loop over trial shape functions
            for j in 1:n_basefuncs
                u  = shape_value(cellvalues, q_point, j)
                ∇u = shape_gradient(cellvalues, q_point, j)
                ## Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
                Me[i, j] += (δu * u) * dΩ
            end
        end
    end
    return nothing
end

function assemble_interface!(Ki::Matrix, iv::InterfaceValues, γ::Real)
    ## Reset to 0
    fill!(Ki, 0)
    ## Loop over quadrature points
    for q_point in 1:getnquadpoints(iv)
        ## Get the normal to face A
        normal = getnormal(iv, q_point)
        ## Get the quadrature weight
        dΓ = getdetJdV(iv, q_point)
        ## Loop over test shape functions
        for i in 1:getnbasefunctions(iv)
            ## Multiply the jump by the normal, as the definition used in Ferrite doesn't include the normals.
            test_jump = shape_value_jump(iv, q_point, i) * normal
            test_grad_avg = shape_gradient_average(iv, q_point, i)

            ## Loop over trial shape functions
            for j in 1:getnbasefunctions(iv)
                ## Multiply the jump by the normal, as the definition used in Ferrite doesn't include the normals.
                trial_jump = shape_value_jump(iv, q_point, j) * normal
                trial_grad_avg = shape_gradient_average(iv, q_point, j)

                ## Add contribution to Ki
                order = Ferrite.getorder(iv.here.fun_values[1].ip)
                dim = Ferrite.getdim(iv.here.fun_values[1].ip)
                Ki[i, j] += -(test_jump ⋅ trial_grad_avg + test_grad_avg ⋅ trial_jump) * dΓ + γ * (test_jump ⋅ trial_jump) * dΓ
            end
        end
    end
    return nothing
end

# function assemble_boundary!(Ke::Matrix, fv::FaceValues, γ::Real)
#     ## Reset to 0
#     fill!(Ke, 0)
#     ## Loop over quadrature points
#     for q_point in 1:getnquadpoints(fv)
#         ## Get the normal to face A
#         normal = getnormal(fv, q_point)
#         ## Get the quadrature weight
#         dΓ = getdetJdV(fv, q_point)
#         ## Loop over test shape functions
#         for i in 1:getnbasefunctions(fv)
#             ## Multiply the jump by the normal, as the definition used in Ferrite doesn't include the normals.
#             test_jump = shape_value(fv, q_point, i) * normal
#             test_grad_avg = shape_gradient(fv, q_point, i)
#             ## Loop over trial shape functions
#             for j in 1:getnbasefunctions(fv)
#                 ## Multiply the jump by the normal, as the definition used in Ferrite doesn't include the normals.
#                 trial_jump = shape_value(fv, q_point, j) * normal
#                 trial_grad_avg = shape_gradient(fv, q_point, j)
#                 order = Ferrite.getorder(fv.fun_values[1].ip)
#                 dim = Ferrite.getdim(fv.fun_values[1].ip)
#                 ## Add contribution to Ki
#                 Ke[i, j] += -(test_jump ⋅ trial_grad_avg + test_grad_avg ⋅ trial_jump) * dΓ + γ * (test_jump ⋅ trial_jump) * dΓ
#             end
#         end
#     end
#     return nothing
# end


function assemble_global(cellvalues::CellValues, interfacevalues::InterfaceValues, M::SparseMatrixCSC, K::SparseMatrixCSC, dh::DofHandler, γ::Real)
    # ∂Ω = union(
    #     getfaceset(grid, "left"),
    #     getfaceset(grid, "right"),
    # );

    ## Allocate the element stiffness matrix and element force vector
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    Me = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    Ki = zeros(n_basefuncs * 2, n_basefuncs * 2)
    ## Allocate global force vector f
    f = zeros(ndofs(dh))
    ## Create an assembler
    assembler_Kf = start_assemble(K, f)
    assembler_M = start_assemble(M)
    ## Loop over all cells
    for cell in CellIterator(dh)
        ## Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        ## Compute volume integral contribution
        assemble_element!(Me, Ke, fe, cellvalues)
        ## Assemble Ke and fe into K and f
        assemble!(assembler_Kf, celldofs(cell), Ke, fe)
        assemble!(assembler_M, celldofs(cell), Me)
    end
    for ic in InterfaceIterator(dh)
        ## Reinitialize interfacevalues for this interface
        reinit!(interfacevalues, ic)
        ## Compute interface surface integrals contribution
        assemble_interface!(Ki, interfacevalues, γ/compute_h(ic))
        ## Assemble Ki into K
        assemble!(assembler_Kf, interfacedofs(ic), Ki)
    end
    # for fc in FaceIterator(dh, ∂Ω)
    #     ## Reinitialize here for this boundary face
    #     reinit!(interfacevalues.here, fc)
    #     ## Compute boundary face surface integrals contribution
    #     assemble_boundary!(Ke, interfacevalues.here, γ)
    #     ## Assemble Ke into K
    #     assemble!(assembler_Kf, celldofs(fc), Ke)
    # end
    return M, K, f
end

interpolate_linear(t₁, t₂, y₁, y₂, t) = t₂>t₁ ? y₁*(t₂-t)/(t₂-t₁) .+ y₂*(t-t₁)/(t₂-t₁) : y₁




struct KellyEstimator
    p::Int
end
getdistance(p1::Vec{N, T},p2::Vec{N, T}) where {N, T} = norm(p1-p2);
#TODO: why this allocates?
getdiameter(cell_coords::AbstractVector{Vec{N, T}}) where {N, T} = maximum(x->(maximum(y->getdistance(y, x), cell_coords; init=0.0)), cell_coords; init=0.0);

Base.@propagate_inbounds function estimate_kelly_interface!(T::Type{<:AbstractFloat}, err::AbstractVector, u::AbstractVector, interface_cache::InterfaceCache, interfacevalues, #=interface_diffusion_cache::BilinearDiffusionInterfaceCache=#)
    error::T = 0.0
    facet_a = interface_cache.a.current_facet_id
    interface_coords = @view getcoordinates(interface_cache)[1][SVector(Ferrite.reference_facets(RefQuadrilateral)[facet_a])]
    # TODO: use getdiameter
    h = getdistance(interface_coords[1], interface_coords[2])
    p = 1
    # @unpack interfacevalues, Dcache_here, Dcache_there = interface_diffusion_cache
    jump_term = 0
    @inbounds for qp in 1:getnquadpoints(interfacevalues.here)
        # D_here = evaluate_coefficient(Dcache_here, interface_cache, qp, time)
        # D_there = evaluate_coefficient(Dcache_there, interface_cache, qp, time)
        normal = getnormal(interfacevalues, qp)
        Wf = getdetJdV(interfacevalues, qp) # maybe?
        # TODO urgent why does jump allocate
        ∇u_jump = function_gradient_jump(interfacevalues, qp, u)
        jump_term += norm2(Wf*(( ∇u_jump ⋅ normal)))
    end
    ret = h/(2p)*jump_term
    err[cellid(interface_cache.a)] += ret
    err[cellid(interface_cache.b)] += ret
end
