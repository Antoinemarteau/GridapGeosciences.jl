@inline _cyl_inv_metric(::Val{D}, r) where D = SymTensorValue{D,Float64,3}(1/r^2, 0.0, 1.0)

struct CylinderInvMetric{D} <: Field
  radius :: Float64
end
CylinderInvMetric(radius::Float64; DG_calc=false) = CylinderInvMetric{DG_calc ? Con(2) : 2}(radius)
Gridap.Arrays.evaluate!(_, m::CylinderInvMetric{D}, x::Point) where D = _cyl_inv_metric(Val(D), m.radius)
function Gridap.Arrays.return_cache(m::CylinderInvMetric{D}, xs::AbstractArray{<:Point}) where D
  CachedArray(fill(_cyl_inv_metric(Val(D), m.radius), size(xs)))
end
function Gridap.Arrays.evaluate!(cache, m::CylinderInvMetric{D}, xs::AbstractArray{<:Point}) where D
  setsize!(cache, size(xs))
  fill!(cache.array, _cyl_inv_metric(Val(D), m.radius))
  cache.array
end

inverse_metric_field(f::CylinderMetric{2})      = CylinderInvMetric{2}(f.radius)
inverse_metric_field(f::CylinderMetric{Cov(2)}) = CylinderInvMetric{Con(2)}(f.radius)
