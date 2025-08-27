function Ferrite.reinit!(iv::InterfaceValues, ic::InterfaceCache{<:FacetCache{<:CellCache{<:Any, <:KoppGrid{Dim}}}}, topology::KoppTopology) where Dim
    root_here = topology.root_idx[cellid(ic.a)]
    root_there = topology.root_idx[cellid(ic.b)]
    facet_here = ic.a.current_facet_id[]
    facet_there = ic.b.current_facet_id[]

    OI_here = topology.root_facet_orientation_info[(root_here-1) *2Dim  + facet_here ]
    OI_there = topology.root_facet_orientation_info[(root_there-1) *2Dim + facet_there ]
    flipped = root_here == root_there ? true : OI_here.flipped ⊻ OI_there.flipped
    shift_index = OI_there.shift_index - OI_here.shift_index
    interface_transformation = Ferrite.InterfaceOrientationInfo{Ferrite.RefHypercube{Dim}, Ferrite.RefHypercube{Dim}}(flipped, shift_index, OI_there.shift_index, facet_here, facet_there)
    return Ferrite.reinit!(iv,
        getcells(ic.a.cc.grid, cellid(ic.a)),
        getcoordinates(ic.a),
        ic.a.current_facet_id[],
        getcells(ic.b.cc.grid, cellid(ic.b)),
        getcoordinates(ic.b),
        ic.b.current_facet_id[],
        interface_transformation
        )
end
function Ferrite.reinit!(
    iv::InterfaceValues,
    cell_here::KoppCell, coords_here::AbstractVector{Vec{dim, T}}, facet_here::Int,
    cell_there::KoppCell, coords_there::AbstractVector{Vec{dim, T}}, facet_there::Int,
    interface_transformation
) where {dim, T}

    # reinit! the here side as normal
    Ferrite.reinit!(iv.here, cell_here, coords_here, facet_here)
    dim == 1 && return reinit!(iv.there, cell_there, coords_there, facet_there)
    # Transform the quadrature points from the here side to the there side
    Ferrite.set_current_facet!(iv.there, facet_there) # Includes boundscheck

    quad_points_a = Ferrite.getpoints(iv.here.fqr, facet_here)
    quad_points_b = Ferrite.getpoints(iv.there.fqr, facet_there)

    temp_a = copy(quad_points_a)
    seq_length = get_refinement_level(cell_here) - get_refinement_level(cell_there)
    facet_a = interface_transformation.facet_a
    flipped = interface_transformation.flipped

    for (idx, point) in pairs(temp_a)
        face_point = Ferrite.element_to_facet_transformation(point, Ferrite.RefHypercube{dim}, facet_a)
        flipped && (face_point *= -1)
        temp_a[idx] = Ferrite.facet_to_element_transformation(face_point, Ferrite.RefHypercube{dim}, facet_a)
    end
    _transform_to_parent!(temp_a, (Ferrite.reference_coordinates(Lagrange{Ferrite.RefHypercube{dim},1}())), cell_here.sequence & ((1 << ((dim + 1) * seq_length)) - 1 ))
    Ferrite.transform_interface_points!(quad_points_b, temp_a, Ferrite.InterfaceOrientationInfo{Ferrite.RefHypercube{dim}, Ferrite.RefHypercube{dim}}(false, interface_transformation.shift_index, interface_transformation.lowest_node_shift_index, interface_transformation.facet_a, interface_transformation.facet_b))
    # TODO: This is the bottleneck, cache it?
    @assert length(quad_points_a) <= length(quad_points_b)

    # Re-evaluate shape functions in the transformed quadrature points
    Ferrite.precompute_values!(Ferrite.get_fun_values(iv.there),  quad_points_b)
    Ferrite.precompute_values!(Ferrite.get_geo_mapping(iv.there), quad_points_b)

    # reinit! the "there" side
    Ferrite.reinit!(iv.there, cell_there, coords_there, facet_there)
    return iv
end
