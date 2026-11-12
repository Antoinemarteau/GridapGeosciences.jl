# ── TorusMap ─────────────────────────────────────────────────────────────────
#
# Chart map of one rectangular patch of a torus, (s,t) ∈ [0,1]² → R³.
# The poloidal angle is θ = 2π(θ₀ + Δθ·g(s)) and the toroidal angle is
# φ = 2π(φ₀ + Δφ·g(t)), with the grading g(s) = s + a·s(1−s) and
# x = ((R + r sinθ) cosφ, (R + r sinθ) sinφ, r cosθ).
#
# g fixes 0 and 1 and has g'(0) = 1+a, g'(1) = 1−a. Two patches sharing an edge
# therefore trace it as the same parametrized curve, while the speed transverse
# to it is (1−a) on one side and (1+a) on the other: the atlas is C⁰ but not C¹
# across patch boundaries unless a = 0.

@inline _torus_grading(a, s) = s + a*s*(1 - s)
@inline _torus_grading_derivative(a, s) = 1 + a*(1 - 2s)

struct TorusMap <: Field
  major_radius :: Float64
  minor_radius :: Float64
  offset       :: Point{2,Float64}   # (θ₀, φ₀), in turns
  width        :: Point{2,Float64}   # (Δθ, Δφ), in turns
  grading      :: Float64            # a, with |a| < 1 for g to be monotone
end

@inline function _torus_angles(m::TorusMap, x)
  a = m.grading
  θ = 2π*(m.offset[1] + m.width[1]*_torus_grading(a, x[1]))
  φ = 2π*(m.offset[2] + m.width[2]*_torus_grading(a, x[2]))
  θ, φ
end

function Gridap.Arrays.evaluate!(_, m::TorusMap, x::Point)
  θ, φ = _torus_angles(m, x)
  ρ = m.major_radius + m.minor_radius*sin(θ)
  Point(ρ*cos(φ), ρ*sin(φ), m.minor_radius*cos(θ))
end

# J[i,j] = ∂xⱼ/∂sᵢ, the layout Gridap uses for the gradient of a vector field.
function Gridap.Arrays.evaluate!(_, f::FieldGradient{1,<:TorusMap}, x::Point)
  m = f.object
  θ, φ = _torus_angles(m, x)
  r, a = m.minor_radius, m.grading
  ρ  = m.major_radius + r*sin(θ)
  θs = 2π*m.width[1]*_torus_grading_derivative(a, x[1])
  φt = 2π*m.width[2]*_torus_grading_derivative(a, x[2])
  TensorValue{2,3,Float64}(
    θs*r*cos(θ)*cos(φ), -φt*ρ*sin(φ),
    θs*r*cos(θ)*sin(φ),  φt*ρ*cos(φ),
    -θs*r*sin(θ),        0.0)
end
