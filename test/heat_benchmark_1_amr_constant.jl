using Ferrite
using Tensors
using Test
using Logging
using ForwardDiff
import SHA
using Random
using LinearAlgebra: norm2, svdvals, pinv, diag
using SparseArrays
using StaticArrays
using OrderedCollections
using WriteVTK
import Metis
using HCubature: hcubature, hquadrature
using Unroll
using DataStructures
import LinearSolve as LS
include("lts_utils.jl")
include("/home/amohamed/.julia/dev/Ferrite/src/Grid/Adaptivity/kopp.jl")
begin
    x = 1.025
    γ = 8.0
    function doassemble_K!(K::SparseMatrixCSC, f::Vector, cellvalues::CellValues, dh::DofHandler)

        n_basefuncs = getnbasefunctions(cellvalues)
        Ke = zeros(n_basefuncs, n_basefuncs)
        fe = zeros(n_basefuncs)

        assembler = start_assemble(K, f)

        for cell in CellIterator(dh)

            fill!(Ke, 0)
            fill!(fe, 0)

            reinit!(cellvalues, cell)

            for q_point in 1:getnquadpoints(cellvalues)
                dΩ = getdetJdV(cellvalues, q_point)

                for i in 1:n_basefuncs
                    v = shape_value(cellvalues, q_point, i)
                    ∇v = shape_gradient(cellvalues, q_point, i)
                    # fe[i] += 0.1 * v * dΩ
                    for j in 1:n_basefuncs
                        ∇u = shape_gradient(cellvalues, q_point, j)
                        Ke[i, j] += 1.0e-2 * (∇v ⋅ ∇u) * dΩ
                    end
                end
            end

            assemble!(assembler, celldofs(cell), Ke, fe)
        end
        return K, f
    end

    function doassemble_M!(M::SparseMatrixCSC, cellvalues::CellValues, dh::DofHandler)

        n_basefuncs = getnbasefunctions(cellvalues)
        Me = zeros(n_basefuncs, n_basefuncs)

        assembler = start_assemble(M)

        for cell in CellIterator(dh)

            fill!(Me, 0)

            reinit!(cellvalues, cell)

            for q_point in 1:getnquadpoints(cellvalues)
                dΩ = getdetJdV(cellvalues, q_point)

                for i in 1:n_basefuncs
                    v = shape_value(cellvalues, q_point, i)
                    for j in 1:n_basefuncs
                        u = shape_value(cellvalues, q_point, j)
                        Me[i, j] += (v * u) * dΩ
                    end
                end
            end

            assemble!(assembler, celldofs(cell), Me)
        end
        return M
    end

    function analytical_solution(x, y, t; D=0.01, a=50)
        denominator = 1 + 4*a*D*t
        return exp(-a*(x^2 + y^2)/denominator) / denominator
    end

    function solve_direct(grid, dh, K, M, u, level, io_enabled::Bool=false) where Dim
        pvd = io_enabled ? paraview_collection("heat-benchmark-1-constant") : nothing
        # dt = 1/maximum(eigen(nMinvK))
        h = 0.1/(2^level)
        dt = 1/3 * (h^2)/(1e-2)
        @info dt
        u_prev = copy(u)
        prob = LS.LinearProblem(M, u_prev)

        linsolve = LS.init(prob)

        for T in 0:0.1:1.9
            io_enabled && VTKGridFile("heat-benchmark-1-constant/heat-$T.vtu", grid; write_discontinuous=true) do vtk
                write_solution(vtk, dh, u, "_")
                pvd[T] = vtk
            end
            t = T
            while t < T + 0.1
                Δt = dt
                if t + Δt > T+ 0.1
                    Δt = T+ 0.1 - t
                end
                Tensors.LinearAlgebra.mul!(u_prev, -K, u)
                linsolve.b .= u_prev
                sol1 = LS.solve!(linsolve)
                u_prev = sol1.u
                u .= u .+ Δt *  u_prev
                # return
                t += Δt
            end
        end
        # notify(uvis)
        io_enabled && vtk_save(pvd)

    end
    refinement_levels = 1:5

    for level in refinement_levels
        refshape = RefQuadrilateral
        grid = generate_grid(Quadrilateral, (10, 10) .* 2^level, Ferrite.Vec((-1.0, -1.0)), Ferrite.Vec((1.0, 1.0)))

        ip = Lagrange{refshape,1}()
        qr = QuadratureRule{refshape}(2)

        face_qr = FacetQuadratureRule{refshape}(2)
        cellvalues = CellValues(qr, ip)
        facetvalues = FacetValues(face_qr, ip)

        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)

        K = allocate_matrix(dh)
        M = allocate_matrix(dh);
        f = zeros(ndofs(dh));
        u = zeros(ndofs(dh))

        K, f = doassemble_K!(K, f, cellvalues, dh)
        M = doassemble_M!(M, cellvalues, dh)
        Ferrite.apply_analytical!(u, dh, :u, x -> analytical_solution(x[1], x[2], 0.0))
        @time solve_direct(grid, dh, K, M, u, level, false)
        @info extrema(u)
    end

end
