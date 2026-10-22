using CairoMakie
using LinearAlgebra: norm

"""
    animate_trail(positions, outpath; kwargs...) -> outpath

Render an MP4 of a single particle moving through a sequence of 2-D positions.
At each frame the particle is drawn as a filled circle and a coloured trail
accumulates behind it, making the trajectory visible throughout the animation.

Arguments
- `positions`: `Vector{Point2f}` — one entry per frame, in order.
- `outpath`:   destination path for the MP4 file (must end in `.mp4`).

Axis limits are inferred from `positions` with `padding` added on each side.
Override either axis with `xlims=(lo, hi)` or `ylims=(lo, hi)` — useful for
enforcing physical constraints such as keeping the ground at y = 0.

Key keyword arguments (with defaults)
- `framerate = 30`             frames per second
- `trail_color = :dodgerblue`  colour of the accumulating path
- `particle_color = :crimson`  colour of the moving dot
- `markersize = 20`            dot diameter in points
- `padding = 2.0`              axis padding when limits are auto-computed
- `xlims / ylims = nothing`    explicit axis bounds; `nothing` means auto
- `figsize = (960, 540)`       figure resolution in pixels
"""
function animate_trail(positions::AbstractVector{Point2f}, outpath::AbstractString;
                       framerate::Int                       = 30,
                       title::AbstractString                = "",
                       xlabel::AbstractString               = "x",
                       ylabel::AbstractString               = "y",
                       trail_color                          = :dodgerblue,
                       particle_color                       = :crimson,
                       markersize::Real                     = 20,
                       padding::Real                        = 2.0,
                       xlims::Union{Nothing,Tuple}          = nothing,
                       ylims::Union{Nothing,Tuple}          = nothing,
                       figsize::Tuple                       = (960, 540))
    xmin = isnothing(xlims) ? minimum(p[1] for p in positions) - padding : xlims[1]
    xmax = isnothing(xlims) ? maximum(p[1] for p in positions) + padding : xlims[2]
    ymin = isnothing(ylims) ? minimum(p[2] for p in positions) - padding : ylims[1]
    ymax = isnothing(ylims) ? maximum(p[2] for p in positions) + padding : ylims[2]

    fig   = Figure(size = figsize)
    ax    = Axis(fig[1, 1]; xlabel, ylabel, title,
                 limits = (xmin, xmax, ymin, ymax))
    trail = Observable(Point2f[])
    pt    = Observable(positions[1])
    lines!(ax, trail; color = trail_color, linewidth = 3)
    scatter!(ax, pt; color = particle_color, markersize)

    record(fig, outpath, eachindex(positions); framerate) do i
        push!(trail[], positions[i]); notify(trail)
        pt[] = positions[i]
    end
    return outpath
end

"""
    plot_function1d(xs, ys; kwargs...) -> (; fig, ax)
    plot_function1d(f, xs; kwargs...) -> (; fig, ax)

Plot a single-variable function y = f(x) as a continuous line.

Two calling forms are supported:
- Pass pre-computed arrays `xs` and `ys` directly.
- Pass a callable `f` and a range `xs`; the function is evaluated as `f.(xs)`.

Returns a named tuple `(; fig, ax)`. Use `ax` to add more curves with `lines!`,
change axis labels, set tick marks, add a legend, etc. Use `fig` to save.

Key keyword arguments (with defaults)
- `color = :dodgerblue`
- `linewidth = 2`
- `figsize = (700, 500)`
"""
function plot_function1d(xs::AbstractVector, ys::AbstractVector;
                         title::AbstractString  = "",
                         xlabel::AbstractString = "x",
                         ylabel::AbstractString = "y",
                         color                  = :dodgerblue,
                         linewidth::Real        = 2,
                         figsize::Tuple         = (700, 500))
    fig = Figure(size = figsize)
    ax  = Axis(fig[1, 1]; xlabel, ylabel, title)
    lines!(ax, xs, ys; color, linewidth)
    return (; fig, ax)
end

plot_function1d(f::Function, xs::AbstractVector; kw...) =
    plot_function1d(xs, f.(xs); kw...)

"""
    plot_polar(θs, rs; kwargs...) -> (; fig, ax)
    plot_polar(f, θs; kwargs...) -> (; fig, ax)

Plot a polar curve r = f(θ) in the Cartesian plane.

Two calling forms are supported:
- Pass angle array `θs` and radius array `rs` directly.
- Pass a callable `f` and an angle range `θs`; evaluated as `f.(θs)`.

The conversion to Cartesian coordinates `(r cos θ, r sin θ)` is handled
internally. The axis uses `DataAspect()` so circles are always round.
Negative radii are allowed and reflect the point through the origin, matching
the standard polar-curve convention.

Returns a named tuple `(; fig, ax)`. Use `ax` to overlay additional curves
(e.g. comparing two polar functions) or to add annotations.

Key keyword arguments (with defaults)
- `color = :dodgerblue`
- `linewidth = 2`
- `figsize = (600, 600)`
"""
function plot_polar(θs::AbstractVector, rs::AbstractVector;
                    title::AbstractString  = "",
                    color                  = :dodgerblue,
                    linewidth::Real        = 2,
                    figsize::Tuple         = (600, 600))
    xs = rs .* cos.(θs)
    ys = rs .* sin.(θs)
    fig = Figure(size = figsize)
    ax  = Axis(fig[1, 1]; title, aspect = DataAspect())
    lines!(ax, xs, ys; color, linewidth)
    return (; fig, ax)
end

plot_polar(f::Function, θs::AbstractVector; kw...) = plot_polar(θs, f.(θs); kw...)

"""
    plot_curve2d(xs, ys; kwargs...) -> (; fig, ax)

Plot a parametric plane curve (x(t), y(t)) as a continuous path.

Unlike `plot_vectorfield2d`, which shows a vector at every grid point,
this function draws only the curve itself — no grid, no arrows. It is the
right choice for trajectories, phase-space orbits, and any path through the
plane that is defined parametrically rather than as y = f(x).

The axis uses `DataAspect()` by default so distances are true to scale
(override with `aspect = nothing` for non-physical axes such as phase space).

Colour the path by a scalar quantity (e.g. speed, time, energy) by passing a
vector of the same length as `xs` to `color_values`. When provided, each line
segment is coloured by the scalar at its start point and a colour bar is added.

Returns a named tuple `(; fig, ax)`. Use `ax` to overlay additional curves
or markers; use `fig` to save or add extra panels.

Key keyword arguments (with defaults)
- `color = :dodgerblue`              used when `color_values` is not given
- `color_values = nothing`           per-point scalar for colour-mapped path
- `colormap = :viridis`              colour map used with `color_values`
- `clabel = ""`                      colour bar label (only with `color_values`)
- `linewidth = 2`
- `aspect = DataAspect()`
- `figsize = (600, 600)`
"""
function plot_curve2d(xs::AbstractVector, ys::AbstractVector;
                      title::AbstractString                       = "",
                      xlabel::AbstractString                      = "x",
                      ylabel::AbstractString                      = "y",
                      color                                       = :dodgerblue,
                      color_values::Union{Nothing,AbstractVector} = nothing,
                      colormap                                    = :viridis,
                      clabel::AbstractString                      = "",
                      linewidth::Real                             = 2,
                      aspect                                      = DataAspect(),
                      figsize::Tuple                              = (600, 600))
    fig = Figure(size = figsize)
    ax  = Axis(fig[1, 1]; xlabel, ylabel, title, aspect)
    if isnothing(color_values)
        lines!(ax, xs, ys; color, linewidth)
    else
        segs = [Point2f(xs[i], ys[i]) => Point2f(xs[i+1], ys[i+1])
                for i in 1:length(xs)-1]
        lc = linesegments!(ax, segs; color = color_values[1:end-1],
                           colormap, linewidth)
        Colorbar(fig[1, 2], lc; label = clabel)
    end
    return (; fig, ax)
end

"""
    plot_heatmap2d(xs, ys, zs; kwargs...) -> (; fig, ax)

Plot a 2-D scalar field as a colour-mapped heat map with an attached colour bar.

`xs` and `ys` are 1-D coordinate vectors; `zs` is a matrix of size
`(length(xs), length(ys))` where `zs[i, j]` is the field value at `(xs[i], ys[j])`.
Typical use cases include temperature distributions, electric potential maps,
probability densities, and any other quantity defined over a 2-D grid.

Returns a named tuple `(; fig, ax)`. Use `ax` to add contour lines, scatter
points, or annotations on top of the heat map; use `fig` to save.

Key keyword arguments (with defaults)
- `colormap = :viridis`  any Makie-compatible colormap name or object
- `clabel = ""`          label on the colour bar
- `figsize = (700, 600)`
"""
function plot_heatmap2d(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix;
                        title::AbstractString  = "",
                        xlabel::AbstractString = "x",
                        ylabel::AbstractString = "y",
                        clabel::AbstractString = "",
                        colormap               = :viridis,
                        figsize::Tuple         = (700, 600))
    fig = Figure(size = figsize)
    ax  = Axis(fig[1, 1]; xlabel, ylabel, title)
    hm  = heatmap!(ax, xs, ys, zs; colormap)
    Colorbar(fig[1, 2], hm; label = clabel)
    return (; fig, ax)
end

"""
    plot_surface3d(xs, ys, zs; kwargs...) -> (; fig, ax)

Plot a 3-D surface z = f(x, y) with an attached colour bar encoding height.

`xs` and `ys` are 1-D coordinate vectors; `zs` is a matrix of size
`(length(xs), length(ys))` where `zs[i, j]` is the height at `(xs[i], ys[j])`.
Rendered with `Axis3`, so the result is a fully 3-D perspective view.
Useful for energy landscapes, wavefunctions, and any height field over a grid.

Returns a named tuple `(; fig, ax)`. `ax` is an `Axis3`; use it to add scatter
points in 3-D space or to adjust the camera elevation and azimuth angle via
`ax.elevation` and `ax.azimuth`.

Key keyword arguments (with defaults)
- `colormap = :viridis`  any Makie-compatible colormap name or object
- `zlabel = "z"`         label for the vertical axis
- `figsize = (800, 650)`
"""
function plot_surface3d(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix;
                        title::AbstractString  = "",
                        xlabel::AbstractString = "x",
                        ylabel::AbstractString = "y",
                        zlabel::AbstractString = "z",
                        colormap               = :viridis,
                        figsize::Tuple         = (800, 650))
    fig = Figure(size = figsize)
    ax  = Axis3(fig[1, 1]; xlabel, ylabel, zlabel, title)
    sf  = surface!(ax, xs, ys, zs; colormap)
    Colorbar(fig[1, 2], sf)
    return (; fig, ax)
end

"""
    plot_vectorfield2d(xs, ys, us, vs; kwargs...) -> (; fig, ax)

Plot a 2-D vector field as a grid of arrows.

`xs` and `ys` are 1-D coordinate vectors defining the grid.
`us` and `vs` are matrices of size `(length(xs), length(ys))` giving the
x- and y-components of the field at each grid point: `us[i, j]` and `vs[i, j]`
correspond to the point `(xs[i], ys[j])`.

By default (`normalize = true`) all arrows are rescaled so the longest arrow
fits within one grid cell, keeping the plot readable regardless of the field's
magnitude range. Set `normalize = false` to preserve the raw magnitudes.
Typical use cases include electric and magnetic fields, fluid flow, and
gradient or curl visualisations.

Returns a named tuple `(; fig, ax)`. Use `ax` to overlay a heat map of the
field magnitude, add streamlines, or mark specific points of interest.

Key keyword arguments (with defaults)
- `arrowcolor = :black`
- `arrowsize = 15`          tip size in points
- `normalize = true`        rescale arrows to fit within one grid cell
- `figsize = (700, 600)`
"""
function plot_vectorfield2d(xs::AbstractVector, ys::AbstractVector,
                            us::AbstractMatrix, vs::AbstractMatrix;
                            title::AbstractString  = "",
                            xlabel::AbstractString = "x",
                            ylabel::AbstractString = "y",
                            arrowcolor             = :black,
                            arrowsize::Real        = 15,
                            normalize::Bool        = true,
                            figsize::Tuple         = (700, 600))
    points = vec([Point2f(x, y) for x in xs, y in ys])
    dirs   = vec([Vec2f(us[i, j], vs[i, j])
                  for i in eachindex(xs), j in eachindex(ys)])

    if normalize && !isempty(dirs)
        dx   = length(xs) > 1 ? abs(xs[2] - xs[1]) : 1.0
        dy   = length(ys) > 1 ? abs(ys[2] - ys[1]) : 1.0
        cell = min(dx, dy)
        mlen = maximum(norm(d) for d in dirs)
        dirs = mlen > 0 ? dirs .* Float32(0.8 * cell / mlen) : dirs
    end

    fig = Figure(size = figsize)
    ax  = Axis(fig[1, 1]; xlabel, ylabel, title)
    arrows!(ax, points, dirs; color = arrowcolor, arrowsize, linewidth = 1.5)
    return (; fig, ax)
end

"""
    savefig(path, fig_or_result; kwargs...) -> path

Save a figure to disk, creating any missing parent directories automatically.

`fig_or_result` accepts either a `Figure` or the named tuple `(; fig, ax)`
returned by any `plot_*` function in this file, so you can pass the result
directly without unpacking first.

The output format is inferred from the file extension (`.png`, `.pdf`, `.svg`).
All keyword arguments are forwarded to CairoMakie's `save`. Useful ones:
- `px_per_unit = 2`  raster resolution multiplier (2 → 2× pixel density for PNG)
- `pt_per_unit = 1`  point-to-pixel ratio for vector formats (PDF, SVG)
"""
function savefig(path::AbstractString, fig::Figure; kw...)
    dir = dirname(path)
    isempty(dir) || mkpath(dir)
    save(path, fig; kw...)
    return nothing
end

savefig(path::AbstractString, result::NamedTuple; kw...) =
    savefig(path, result.fig; kw...)