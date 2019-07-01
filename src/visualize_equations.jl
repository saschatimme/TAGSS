"""
    compute_z(f, rx, ry)

Compute all values ``f(x,y)`` for ``(x,y) ∈ rx × ry``.
"""
compute_z(f::DP.Polynomial, rx, ry) = compute_z(SP.Polynomial(f), rx, ry)
compute_z(f::SP.Polynomial, rx, ry) = [f(SVector(x,y)) for x in rx, y in ry]
compute_z(f::Function, rx, ry) = [f(x,y) for x in rx, y in ry]

"""
    curvature(f, xs, ys)

Compute the curvature of `f` for the diagonal of `xs` and `ys`.
"""
function curvature(f::SP.Polynomial, xs::AbstractVector, ys::AbstractVector)
    map((x,y) -> curvature(f, SVector(x,y)), xs, ys)
end
function curvature(f::SP.Polynomial, x)
    ∇f = SP.gradient(f, x)
    H = SP.hessian(f, x)
    P = I - (∇f * ∇f') / norm(∇f)^2
    Q = qr(P).Q[:,1:end-1]
    E = eigvals(Symmetric(Q'*H*Q))
    _, idx = findmax(abs.(E))
    -E[idx] / norm(∇f)
end



function visualize(f::DP.Polynomial; scene_resolution=(1000,1000), kwargs...)
    scene = Scene(resolution=scene_resolution, scale_plot=false)
    visualize!(scene, f; kwargs...)
end
function visualize!(scene, f::DP.Polynomial; wireframe=nothing, kwargs...)
    if DP.nvariables(f) == 2
        implicit_curve!(scene, f; kwargs...)
    else
        implicit_surface!(scene, f; wireframe=wireframe, kwargs...)
    end
end

"""
    implicit_curve!(scene, f; x_min=-5, x_max=5, y_min=x_min, y_max=x_max,
                      color_curvature=false, color=:steelblue, resolution=1000)

Visualize the implicit curve `f` in the box `[x_min, x_max] × [y_min, y_max]`.
If `color_curvature` is `true` then the curve is locally colored depending on its curvature.
Otherwise `color` is used.
"""
function implicit_curve!(scene, f::DP.Polynomial; x_min=-3., x_max=3., y_min=x_min, y_max=x_max,
            color_curvature=false, color=:steelblue, resolution=1000, grid=true, kwargs...)
    g = SP.Polynomial(f)
    x = collect(range(x_min, x_max; length=resolution))
    y = collect(range(y_min, y_max; length=resolution))
    z = compute_z(f, x, y)

    limits = FRect(x_min, y_min, x_max - x_min, y_max - y_min)
    lvl = Contour.contour(x, y, z, 0.0)
    lines = Contour.lines(lvl)
    !isempty(lines) || return scene
    for l in lines
        xs, ys = Contour.coordinates(l)
        if color_curvature
            t = log.(abs.(curvature(g, xs, ys)))
            lines!(scene, xs, ys; linewidth=3, color=t, limits=limits, kwargs...)
        else
            lines!(scene, xs, ys; linewidth=3, color=color, limits=limits, kwargs...)
        end
    end

    axis = scene[Axis] # get the axis object from the scene
    if !isnothing(axis)
        if grid == false
            axis[:grid][:linewidth] = (0,0)
        end
        axis[:ticks][:textsize] = (2, 2)
        axis[:ticks][:textsize] = (2, 2)
        axis[:names][:textsize] = (2, 2)
    end
    scene
end

"""
    implicit_surface!(scene, f; x_min=-3, x_max=3, y_min=x_min, y_max=x_max, z_min=x_min, z_max=z_max,
                      color_curvature=false, color=:steelblue, resolution=1000)

Visualize the implicit curve `f` in the box `[x_min, x_max] × [y_min, y_max]`.
If `color_curvature` is `true` then the curve is locally colored depending on its curvature.
Otherwise `color` is used.
"""
function implicit_surface!(scene, f::DP.Polynomial; x_min=-3., x_max=3., y_min=x_min, y_max=x_max,
                z_min=x_min, z_max=x_max, color_curvature=false, color=:steelblue,
                show_axis=true,
                wireframe=false,
                mesh_resolution=0.04, grid=true, kwargs...)
    g = SP.Polynomial(f)
    box = HyperRectangle(Vec(x_min,y_min, float(z_min)),Vec(float(x_max-x_min),y_max-y_min,z_max-z_min))
    sdf = SignedDistanceField(v -> g(v), box, mesh_resolution)
    m = GLNormalMesh(sdf, MarchingTetrahedra())
    curvature_colors = log.(abs.(curvature.(Ref(g), m.vertices)))
    if isnothing(wireframe) || !wireframe
        mesh!(scene, m, show_axis=show_axis,
                    color=color_curvature ? curvature_colors : color,
                    shading=!color_curvature,
                colormap=ColorSchemes.YlGnBu_5.colors)
    else
        wireframe!(scene, m, show_axis=show_axis, linewidth=0.02, transparency=true, alpha=0.1)
    end

    scene
end
