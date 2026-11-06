function set_num_points!(data::PointVortexData, Np)
    resize_no_copy!(data.points, Np)
    resize_no_copy!(data.signs, Np)
    data
end

function add_point_charges!(data::PointVortexData{T}, vortices::PointVortexState{T}, params::ParamsPointVortex{T}) where {T}
    (; Ls) = params
    (; points, signs) = vortices
    Base.require_one_based_indexing(points)
    Np = length(points)
    set_num_points!(data, Np)
    Threads.@threads for i in eachindex(points, signs)
        @inbounds data.points[i] = Filaments.fold_coordinates_periodic(points[i], Ls)
        @inbounds data.signs[i] = signs[i]
    end
    data
end
