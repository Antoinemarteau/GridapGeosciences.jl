# ── CylinderMetric ────────────────────────────────────────────────────────────
#
# Pullback metric g = diag(r², 1) — constant in chart coordinates.
# Inverse g⁻¹ = diag(1/r², 1).

@inline _cyl_metric(::Val{D}, r) where D = SymTensorValue{D,Float64,3}(r^2, 0.0, 1.0)

struct CylinderMetric{D} <: Field
  radius :: Float64
end
CylinderMetric(radius::Float64; DG_calc=false) = CylinderMetric{DG_calc ? Cov(2) : 2}(radius)
Gridap.Arrays.evaluate!(_, m::CylinderMetric{D}, x::Point) where D = _cyl_metric(Val(D), m.radius)
function Gridap.Arrays.return_cache(m::CylinderMetric{D}, xs::AbstractArray{<:Point}) where D
  CachedArray(fill(_cyl_metric(Val(D), m.radius), size(xs)))
end
function Gridap.Arrays.evaluate!(cache, m::CylinderMetric{D}, xs::AbstractArray{<:Point}) where D
  setsize!(cache, size(xs))
  fill!(cache.array, _cyl_metric(Val(D), m.radius))
  cache.array
end
