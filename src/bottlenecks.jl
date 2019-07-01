export bottlenecks, bottleneck_equations, visualize_bottlenecks, bottleneck_animation

"""
    bottleneck_equations(f)

Create the bottleneck equations. This assumes ``f`` is a complete intersection.
"""
bottleneck_equations(f::DP.AbstractPolynomial) = bottleneck_equations([f])
function bottleneck_equations(f::Vector{<:DP.AbstractPolynomial})
    x = DP.variables(f)
    n, m = length(x), length(f)
    DP.@polyvar y[1:n] v[1:m] w[1:m]
    J = [DP.differentiate(fᵢ, xᵢ) for fᵢ in f, xᵢ in x]
    f′ = [DP.subs(fᵢ, x=>y) for fᵢ in f]
    J′ = [DP.subs(gᵢ, x=>y) for gᵢ in J]
    Nx = (x - y) - J' * v
    Ny = (x - y) - J′' * w
    I = [f; f′; Nx; Ny]
    I
end


"""
    bottlenecks(f)

Compute all bottlenecks of ``f``.
"""
function bottlenecks(result::HC.Result, dim::Integer)
    real_sols = HC.real_solutions(result; only_nonsingular=true)
    duplicate_coordinates = map(s -> s[1:2*dim], real_sols)
    relabel = s -> [s[dim+1:2dim]; s[1:dim]]
    points = HC.unique_points(duplicate_coordinates; group_action=relabel)
    map(p -> (p[1:dim], p[dim+1:2dim]), points)
end
bottlenecks(f::DP.AbstractPolynomial; kwargs...) = bottlenecks([f]; kwargs...)
function bottlenecks(f::Vector{<:DP.AbstractPolynomial}; kwargs...)
    n = DP.nvariables(f)
    res = HC.solve(bottleneck_equations(f); kwargs...)
    bottlenecks(res, n)
end

"""
    visualize_bottlenecks(f, bn_pts)
"""
function visualize_bottlenecks(f, bn_pts::Vector{<:NTuple{2}}; only_minimal=false, kwargs...)
    scene = visualize(f; wireframe=true, kwargs...)
    visualize_bottlenecks!(scene, bn_pts)
end
function visualize_bottlenecks!(scene, bn_pts::Vector{<:NTuple{2}}; only_minimal=false)
    !isempty(bn_pts) || return scene
    n = length(bn_pts[1][1])
    bn_distance = map(((p,q),) -> norm(p-q), bn_pts)
    bn_distance_sort_ind = sortperm(bn_distance)
    n_min_bn = let
        d = bn_distance[bn_distance_sort_ind[1]]
        k = 1
        while k + 1 ≤ length(bn_pts) &&
              bn_distance[bn_distance_sort_ind[k+1]] - d ≤ 1e-8
            k += 1
        end
        k
    end
    # draw surface
    w = maximum(scene.limits.val.widths)
    bn_segments = map(n_min_bn+1:length(bn_pts)) do i
        to_point_pair(bn_pts[bn_distance_sort_ind[i]]...)
    end
    if !isempty(bn_segments) && !only_minimal
        linesegments!(scene, bn_segments; color=:SLATEGRAY, linewidth=2,
                    linestyle=n == 2 ? :dash : nothing, show_axis=false)
    end

    # index of minimal bottleneck
    for i in 1:n_min_bn
        p, q = bn_pts[bn_distance_sort_ind[i]]
        # scatter!(scene, [p[1], q[1]], [p[2], q[2]], markersize=w*0.01, color=:INDIANRED, show_axis=false)
        linesegments!(scene, [to_point_pair(p,q)], linewidth=7.5, color=:INDIANRED, show_axis=false)
    end
    scene
end

function bottleneck_animation(f, bn_pts, filename; resolution=(1000,1000), show_axis=false, seconds=12, kwargs...)
    scene = Scene(resolution=resolution);
    visualize!(scene, f; show_axis=show_axis, color_curvature=true, kwargs...)
    visualize!(scene, f; show_axis=show_axis, wireframe=true, kwargs...)
    visualize_bottlenecks!(scene, bn_pts)

    push!(scene.plots[1].attributes.visible, true)
    push!(scene.plots[2].attributes.visible, false)
    push!(scene.plots[3].attributes.visible, false)
    push!(scene.plots[4].attributes.visible, false)
    N = seconds * 24
    record(scene, filename, 1:N) do i
        if i == div(N,4)
            push!(scene.plots[1].attributes.visible, false)
            push!(scene.plots[2].attributes.visible, true)
        elseif i == 2div(N,4)
            push!(scene.plots[3].attributes.visible, true)
            push!(scene.plots[4].attributes.visible, true)
        end

        rotate_cam!(scene, 4π/N, 0.00, 0.00)
    end
end
