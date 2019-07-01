export ed_equations, ed, solve_ed, visualize_ed, visualize_ed!

"""
    ed_equations(f)

Create the equations to compute all critical points of the euclidean distance function
between a point `u` and V(f).
Returns the equations and the variables corresponding to the point `u` for which
the problem should be computed.
This assumes that `f` is a complete intersection.
"""
ed_equations(f::DP.AbstractPolynomial) = ed_equations([f])
function ed_equations(f::Vector{<:DP.AbstractPolynomial})
    x = DP.variables(f)
    n, m = length(x), length(f)
    DP.@polyvar λ[1:m] u[1:n]
    J = [DP.differentiate(fᵢ, xᵢ) for fᵢ in f, xᵢ in x]
    Nx = (x - u) - J' * λ # lagrange multiplier
    [f; Nx], u
end

"""
    ed(f, u₀)

Compute all critical points of the euclidean distance function
between `u₀` and V(f).
"""
function ed(result::HC.Result, dim::Integer)
    realsols = HC.real_solutions(result; only_nonsingular=true)
    map(p -> p[1:dim], realsols)
end
ed(f::DP.AbstractPolynomial, u₀; kwargs...) = ed([f], u₀; kwargs...)
function ed(f::Vector{<:DP.AbstractPolynomial}, u₀::Vector{<:Number}; kwargs...)
    n = DP.nvariables(f)
    ed(solve_ed(f, u₀), n)
end

solve_ed(f::DP.AbstractPolynomial, u₀; kwargs...) = solve_ed(f, u₀; kwargs...)
function solve_ed(f::Vector{<:DP.AbstractPolynomial}, u₀::Vector{<:Number}; kwargs...)
    F, u = ed_equations(f)
    F₀ = DP.subs.(F, Ref(u => u₀))
    HC.solve(F₀; kwargs...)
end

"""
    visualize_ed(f, u₀, ed_pts)

Visualize all computed solutions of the ed problem. `f` has to be a plane curve
or a surface in three space.
"""
visualize_ed(f, u₀; kwargs...) = visualize_ed(f, u₀, ed(f, u₀); kwargs...)
function visualize_ed(f, u₀, ed_pts::Vector; show_axis=true, kwargs...)
    scene = visualize(f; wireframe=true, show_axis=show_axis, kwargs...)
    visualize_ed!(scene, u₀, ed_pts; show_axis=show_axis)
end
function visualize_ed!(scene, u₀, ed_pts::Vector; show_axis=true)
    if length(u₀) == 2
        scatter!(scene, u₀[1:1], u₀[2:2], markersize=0.05, color=:INDIANRED, show_axis=show_axis)
    else
        scatter!(scene, u₀[1:1], u₀[2:2], u₀[3:3], markersize=0.05, color=:INDIANRED, show_axis=show_axis)
    end

    !isempty(ed_pts) || return scene

    ed_distance = map(p -> norm(p-u₀), ed_pts)

    # find all points of minimal distance, this is a little bit more involved
    # to handle the case that there a multiple minimal points
    # otherwise you could just do
    # val, idx = findmin(ed_distance)
    sorted_ed_pts = ed_pts[sortperm(ed_distance)]
    n_min_ed = let
        k = 1
        while k + 1 ≤ length(ed_pts) &&
            norm(sorted_ed_pts[k+1] - sorted_ed_pts[1]) ≤ 1e-8
            k += 1
        end
        k
    end

    ed_segments = map(p -> to_point_pair(p,u₀), sorted_ed_pts[n_min_ed+1:end])
    linesegments!(scene, ed_segments; color=:INDIANRED, linewidth=2,
        linestyle=:dash, show_axis=show_axis)

    min_ed_segments = map(p -> to_point_pair(p,u₀), sorted_ed_pts[1:n_min_ed])
    linesegments!(scene, min_ed_segments; color=:INDIANRED, linewidth=2,
                show_axis=show_axis)
    if length(u₀) == 2
        scatter!(scene, Point2f0.(ed_pts), markersize=0.05, marker=:x, color=:INDIANRED, show_axis=show_axis)
    else
        scatter!(scene, Point3f0.(ed_pts), markersize=0.05, marker=:x, color=:INDIANRED, show_axis=show_axis)
    end
    scene
end

function to_point_pair(p, q)
    if length(p) == 2
        Point2f0(p) => Point2f0(q)
    else
        Point3f0(p) => Point3f0(q)
    end
end
