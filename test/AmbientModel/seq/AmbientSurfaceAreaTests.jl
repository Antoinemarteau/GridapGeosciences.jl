include("../AmbientSurfaceArea.jl")

## Serial model: 2D
n_ref_lvls = 4
radii = [1.0, 2.0]

for radius in radii
  ambient_models = generate_refined_models(n_ref_lvls, CubedSphereMesh(radius), ExtrinsicManifold())
  parametric_models = generate_refined_models(n_ref_lvls, CubedSphereMesh(radius), IntrinsicManifold())
  AmbientSurfaceArea.main(parametric_models, ambient_models)
end

## Volume form: the parametric metric is read by `⋆` only with DG_calc=true
for radius in radii
  ambient_models = generate_refined_models(n_ref_lvls, CubedSphereMesh(radius), ExtrinsicManifold())
  parametric_models = generate_refined_models(n_ref_lvls, CubedSphereMesh(radius),
                                              IntrinsicManifold(); DG_calc=true)
  AmbientSurfaceArea.main_volume_form(parametric_models, ambient_models)
end

