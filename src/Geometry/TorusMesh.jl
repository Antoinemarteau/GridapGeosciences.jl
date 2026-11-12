# ============================================================
# TorusMesh
# ============================================================

"""
    TorusMesh(major_radius=5.0, minor_radius=1.0; n_poloidal=3, n_toroidal=6, grading=0.0)

`n_poloidal × n_toroidal` atlas of a torus. Every chart is the unit square
`[0,1]²`, mapped by `TorusMap` onto one rectangle of the (θ, φ) angle torus.

With `grading = a ≠ 0` each chart is reparametrized by `g(s) = s + a·s(1−s)` in
both directions. Neighbouring charts still agree on their shared edge as a
parametrized curve, so refined nodes on it coincide, but the parametrization is
only C⁰ across it: the transverse speed is `1−a` on one side and `1+a` on the
other. With `a = 0` the atlas is a uniform grid in (θ, φ), C^∞ across charts.

Each periodic direction needs at least 3 charts: Gridap identifies edges by
their unordered vertex pair, and with 2 cells in a periodic row both cells
would share both of their edges in that direction (see `CylinderMesh`).

The metric is the pullback `JᵀJ` of the ambient map, so `DG_calc=true` is not
supported.
"""
struct TorusMesh <: CoarseMesh
  major_radius :: Float64
  minor_radius :: Float64
  n_poloidal   :: Int
  n_toroidal   :: Int
  grading      :: Float64
  function TorusMesh(major_radius=5.0, minor_radius=1.0;
                     n_poloidal=3, n_toroidal=6, grading=0.0)
    n_poloidal ≥ 3 && n_toroidal ≥ 3 ||
      throw(ArgumentError("TorusMesh needs at least 3 charts per direction, got $n_poloidal × $n_toroidal"))
    abs(grading) < 1 ||
      throw(ArgumentError("TorusMesh grading must satisfy |a| < 1 for the charts to be invertible, got $grading"))
    new(major_radius, minor_radius, n_poloidal, n_toroidal, grading)
  end
end

"""
    get_coarse_mesh(m::TorusMesh) → CoarseMeshInfo{2}

Periodic `n_poloidal × n_toroidal` QUAD mesh. Node `(i,j)`, `0 ≤ i < n_poloidal`,
`0 ≤ j < n_toroidal`, has id `1 + i + n_poloidal*j`; cell `(i,j)` has corners
`(i,j), (i+1,j), (i,j+1), (i+1,j+1)` taken periodically, in Gridap's Z-order
BL, BR, TL, TR, so that every chart's `s` runs along θ and `t` along φ.
"""
function get_coarse_mesh(m::TorusMesh; DG_calc=false)
  DG_calc && throw(ArgumentError("TorusMesh only provides the pullback metric; DG_calc=true is not supported"))
  nθ, nφ = m.n_poloidal, m.n_toroidal
  node(i, j) = Int32(1 + mod(i, nθ) + nθ*mod(j, nφ))

  # Coordinate values are unused — AtlasGrid replaces them with cell_chart_coords.
  node_coords = [Point(Float64(i), Float64(j)) for j in 0:nφ-1 for i in 0:nθ-1]

  cells = [(i, j) for j in 0:nφ-1 for i in 0:nθ-1]
  cell_node_data = Int32[n for (i, j) in cells
                           for n in (node(i,j), node(i+1,j), node(i,j+1), node(i+1,j+1))]
  cell_node_ptrs = Int32.(1:4:4*length(cells)+1)
  cell_node_ids  = Gridap.Arrays.Table(cell_node_data, cell_node_ptrs)
  cell_types     = fill(Int32(1), length(cells))
  reffe          = Gridap.ReferenceFEs.LagrangianRefFE(Float64, QUAD, 1)
  grid = Gridap.Geometry.UnstructuredGrid(
    node_coords, cell_node_ids, [reffe], cell_types, Gridap.Geometry.Oriented())

  topo   = Gridap.Geometry.UnstructuredGridTopology(grid)
  labels = Gridap.Geometry.FaceLabeling(topo)   # closed manifold: no boundary
  model  = Gridap.Geometry.UnstructuredDiscreteModel(grid, topo, labels)

  unit_square = [Point(0.0,0.0), Point(1.0,0.0), Point(0.0,1.0), Point(1.0,1.0)]
  cell_chart_coords = fill(unit_square, length(cells))
  width = Point(1/nθ, 1/nφ)
  ambient_maps = [TorusMap(m.major_radius, m.minor_radius, Point(i/nθ, j/nφ), width, m.grading)
                  for (i, j) in cells]
  metric_fields = _pullback_metrics(ambient_maps)

  CoarseMeshInfo(model, cell_chart_coords, ambient_maps, metric_fields)
end
