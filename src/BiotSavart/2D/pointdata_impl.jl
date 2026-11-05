function set_num_points!(data::PointVortexData, Np)
    resize_no_copy!(data.points, Np)
    resize_no_copy!(data.charges, Np)
    data
end

function add_point_charges!(data::PointVortexData, xp::PointVortices, charges::AbstractVector, params::ParamsPointVortex)
    Base.require_one_based_indexing(xp)
    Np = length(xp)
    set_num_points!(data, Np)
    Threads.@threads for i in eachindex(xp, charges)
        @inbounds data.points[i] = xp[i]
        @inbounds data.charges[i] = charges[i]
    end
    data
end
