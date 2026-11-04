@inline function _mobius_inv_metric(::Val{D}, R, W, offset, x) where D
  θ   = π*(x[1] + offset)/2
  ρ   = R + W*x[2]*cos(θ/2)
  g11 = (π/4)^2 * W^2 * x[2]^2 + (π/2)^2 * ρ^2
  SymTensorValue{D,Float64,3}(1/g11, 0.0, 1/W^2)
end

struct MobiusInvMetric{D} <: Field
  radius :: Float64; half_width :: Float64; theta_offset :: Float64
end
MobiusInvMetric(R::Float64, W::Float64, offset::Float64; DG_calc=false) =
  MobiusInvMetric{DG_calc ? Con(2) : 2}(R, W, offset)
Gridap.Arrays.evaluate!(_, m::MobiusInvMetric{D}, x::Point) where D =
  _mobius_inv_metric(Val(D), m.radius, m.half_width, m.theta_offset, x)
function Gridap.Arrays.return_cache(_::MobiusInvMetric{D}, xs::AbstractArray{<:Point}) where D
  CachedArray(similar(xs, SymTensorValue{D,Float64,3}))
end
function Gridap.Arrays.evaluate!(cache, m::MobiusInvMetric{D}, xs::AbstractArray{<:Point}) where D
  setsize!(cache, size(xs))
  R, W, offset = m.radius, m.half_width, m.theta_offset
  @inbounds for i in eachindex(xs)
    cache.array[i] = _mobius_inv_metric(Val(D), R, W, offset, xs[i])
  end
  cache.array
end

inverse_metric_field(f::MobiusMetric{2}) =
  MobiusInvMetric{2}(f.radius, f.half_width, f.theta_offset)
inverse_metric_field(f::MobiusMetric{Cov(2)}) =
  MobiusInvMetric{Con(2)}(f.radius, f.half_width, f.theta_offset)
