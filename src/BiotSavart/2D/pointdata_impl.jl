function set_num_points!(data::PointVortexData, Np)
    resize_no_copy!(data.points, Np)
    resize_no_copy!(data.signs, Np)
    data
end

function add_point_charges!(data::PointVortexData, xp::PointVortices, signs::AbstractVector, params::ParamsPointVortex)
    (; Ls) = params
    Base.require_one_based_indexing(xp)
    Np = length(xp)
    set_num_points!(data, Np)
    Threads.@threads for i in eachindex(xp, signs)
        @inbounds data.points[i] = Filaments.fold_coordinates_periodic(xp[i], Ls)
        @inbounds data.signs[i] = signs[i]
    end
    data
end
