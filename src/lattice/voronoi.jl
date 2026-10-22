# Voronoi-cell / polyhedron utilities for Wigner-Seitz and Brillouin-zone
# construction.  Two-dimensional cells use Sutherland-Hodgman polygon
# clipping; three-dimensional cells use vertex enumeration from triples
# of bounding planes.

using LinearAlgebra

# ── 2-D ───────────────────────────────────────────────────────────────────────

"""
    clip_polygon(poly, normal, offset) → poly_clipped

Sutherland-Hodgman: keep vertices satisfying  normal · x ≤ offset.
`poly` is a `Vector{NTuple{2,Float64}}` in CCW order.
"""
function clip_polygon(poly::Vector{NTuple{2,Float64}},
                      normal::NTuple{2,Float64}, offset::Real;
                      tol::Real=1e-12)
    n = length(poly);  n == 0 && return poly
    out = NTuple{2,Float64}[]
    for i in 1:n
        p1 = poly[i]
        p2 = poly[mod1(i+1, n)]
        d1 = normal[1]*p1[1] + normal[2]*p1[2] - offset
        d2 = normal[1]*p2[1] + normal[2]*p2[2] - offset
        in1 = d1 ≤ tol;  in2 = d2 ≤ tol
        if in1 && in2
            push!(out, p2)
        elseif in1 && !in2
            t = d1 / (d1 - d2)
            push!(out, (p1[1] + t*(p2[1]-p1[1]), p1[2] + t*(p2[2]-p1[2])))
        elseif !in1 && in2
            t = d1 / (d1 - d2)
            push!(out, (p1[1] + t*(p2[1]-p1[1]), p1[2] + t*(p2[2]-p1[2])))
            push!(out, p2)
        end
    end
    out
end

"""
    voronoi_cell_2d(p, neighbors; bbox=(-10,10,-10,10)) → Vector{NTuple{2,Float64}}

Voronoi cell of point `p` against `neighbors` (all other relevant sites),
returned as the CCW vertex list of the bounded polygon.  Bisector with
neighbor `q` enforces  (q-p) · x ≤ (|q|² − |p|²)/2.
"""
function voronoi_cell_2d(p::NTuple{2,Float64},
                          neighbors::AbstractVector{<:NTuple{2,Real}};
                          bbox::NTuple{4,Real}=(-10.0,10.0,-10.0,10.0))
    xmin, xmax, ymin, ymax = bbox
    poly = NTuple{2,Float64}[(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)]
    for q in neighbors
        (q[1] == p[1] && q[2] == p[2]) && continue
        nx = q[1] - p[1];  ny = q[2] - p[2]
        offset = ((q[1]^2 + q[2]^2) - (p[1]^2 + p[2]^2)) / 2
        poly = clip_polygon(poly, (nx, ny), offset)
        isempty(poly) && return poly
    end
    poly
end

"""
    ws_cell_2d(neighbors; bbox) → Vector{NTuple{2,Float64}}

Wigner-Seitz cell at the origin: the Voronoi cell of (0,0) given neighbor
lattice points.  Bisector becomes  R · x ≤ |R|²/2.
"""
function ws_cell_2d(neighbors::AbstractVector{<:NTuple{2,Real}};
                     bbox::NTuple{4,Real}=(-10.0,10.0,-10.0,10.0))
    voronoi_cell_2d((0.0, 0.0), neighbors; bbox)
end

# ── 3-D ───────────────────────────────────────────────────────────────────────

"""
    ws_cell_3d(neighbors; tol=1e-7)
        → (vertices::Vector{Vector{Float64}}, faces::Vector{Vector{Int}})

Wigner-Seitz cell at the origin in 3-D.  Each `neighbor` R defines the half-
space  R · x ≤ |R|²/2.  Algorithm:

  1. Enumerate triples (i, j, k) of half-spaces; solve the 3×3 linear
     system to get the candidate vertex v.
  2. Discard v if any other half-space is violated (v not in the polytope).
  3. Deduplicate within `tol`.
  4. Each vertex lies on one or more bounding planes; collect, for each
     plane f, the set of vertices on f and sort them CCW around the face
     centroid in the plane's local 2-D frame.

Returns the vertex list and a `faces` vector of per-face CCW vertex indices.
Empty faces (< 3 vertices) are kept so face indices align with `neighbors`.
"""
function ws_cell_3d(neighbors::AbstractVector;
                     tol::Real=1e-7)
    N = length(neighbors)
    normals = [Float64.(collect(R)) for R in neighbors]
    offsets = [dot(n, n) / 2 for n in normals]

    vertices       = Vector{Vector{Float64}}()
    vertex_planes  = Vector{Vector{Int}}()

    @inbounds for i in 1:N, j in (i+1):N, k in (j+1):N
        A = vcat(normals[i]', normals[j]', normals[k]')
        detA = det(A)
        abs(detA) < 1e-10 && continue
        v = A \ [offsets[i], offsets[j], offsets[k]]

        # Reject if any half-space violated
        keep = true
        for m in 1:N
            if dot(normals[m], v) > offsets[m] + tol
                keep = false;  break
            end
        end
        keep || continue

        # Find all planes this vertex lies on
        on_planes = Int[]
        for m in 1:N
            abs(dot(normals[m], v) - offsets[m]) < tol && push!(on_planes, m)
        end

        # Deduplicate
        dup = false
        for w in vertices
            if norm(v .- w) < tol
                dup = true;  break
            end
        end
        dup && continue

        push!(vertices, v)
        push!(vertex_planes, on_planes)
    end

    # Group vertices by face
    face_verts = [Int[] for _ in 1:N]
    for (vi, planes) in enumerate(vertex_planes)
        for f in planes
            push!(face_verts[f], vi)
        end
    end

    # Sort each face CCW about its centroid
    faces = Vector{Vector{Int}}(undef, N)
    for f in 1:N
        idx = face_verts[f]
        if length(idx) < 3
            faces[f] = Int[]
            continue
        end
        n̂ = normals[f] / norm(normals[f])
        centroid = sum(vertices[i] for i in idx) / length(idx)
        # in-plane basis (e1, e2)
        tmp = abs(n̂[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
        e1  = cross(n̂, tmp);  e1 ./= norm(e1)
        e2  = cross(n̂, e1)
        angles = [atan(dot(vertices[i] .- centroid, e2),
                        dot(vertices[i] .- centroid, e1))  for i in idx]
        faces[f] = idx[sortperm(angles)]
    end

    vertices, faces
end

"""
    edges_from_faces(faces) → Vector{Tuple{Int,Int}}

Unique undirected edges, derived from face vertex loops.  Useful for
plotting a wireframe polyhedron with `linesegments!`.
"""
function edges_from_faces(faces::AbstractVector)
    s = Set{Tuple{Int,Int}}()
    for f in faces
        n = length(f);  n < 2 && continue
        for i in 1:n
            a = f[i];  b = f[mod1(i+1, n)]
            push!(s, (min(a,b), max(a,b)))
        end
    end
    collect(s)
end
