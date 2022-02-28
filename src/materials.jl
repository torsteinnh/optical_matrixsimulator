module materials

using Interpolations
using DelimitedFiles


export LoadMaterial, spesifics


function LoadMaterial(path) # Loads material where complex refractive index is contained in one file
    material = readdlm(path, ',', header=true)[1][:, 1:3]

    n_data = [material[:, 1] .* 1e-6, material[:, 2]]
    k_data = [material[:, 1] .* 1e-6, material[:, 3]]
    n_real = LinearInterpolation(n_data[1], n_data[2])
    n_complex = LinearInterpolation(k_data[1], k_data[2])

    n_total(λ) = n_real(λ) - abs(n_complex(λ)) * 1im
    n_total
end

function LoadMaterial(path_n, path_k, number) # Loads material where n and k are split and dumped by webplotdigitizer
    material_n = readdlm(path_n, ',', skipstart=2)[:, (number*2 - 1):(number*2)]
    material_k = readdlm(path_k, ',', skipstart=2)[:, (number*2 - 1):(number*2)]

    material_n = material_n[(material_n[:, 1] .!= ""), :]
    material_k = material_k[(material_k[:, 1] .!= ""), :]

    material_n = material_n[sortperm(material_n[:, 1]), :]
    material_k = material_k[sortperm(material_k[:, 1]), :]

    n_estimator = LinearInterpolation(material_n[:, 1] .* 1e-9, material_n[:, 2])
    k_estimator = LinearInterpolation(material_k[:, 1] .* 1e-9, material_k[:, 2])

    n_total(λ) = n_estimator(λ) - abs(k_estimator(λ)) * 1im
    n_total
end


module spesifics

import ..LoadMaterial
using ...analyticalmaterials


export Air, Au_Johnson, Au_Werner, Au_11nm_Rosenblatt, Au_21nm_Rosenblatt, Au_44nm_Rosenblatt, Au_Babar, Au_McPeak, Au_OlmonEvaporated, Au_OlmonSingleChrystaline, Au_OlmonTemplateStripped, Pd_Johnson, Pd_Werner, Pd_Palm_2018, Ag, LiF, SiO2_core_Sellmeier, SiO2_thinfilm_Ciprian, Au_unloaded, Pd014_unloaded, Pd034_unloaded, Pd034_loaded, Pd042_unloaded, Pd042_loaded, Pd052_unloaded, Pd052_loaded, Pd073_unloaded, Pd073_loaded, Pd_unloaded, Pd_loaded

Air(λ) = 1


Au_11nm_Rosenblatt_estimator = LoadMaterial("materials/refractive_index/Au_11nm_Rosenblatt.csv")
Au_11nm_Rosenblatt(λ) = Au_11nm_Rosenblatt_estimator(λ)

Au_21nm_Rosenblatt_estimator = LoadMaterial("materials/refractive_index/Au_21nm_Rosenblatt.csv")
Au_21nm_Rosenblatt(λ) = Au_21nm_Rosenblatt_estimator(λ)

Au_44nm_Rosenblatt_estimator = LoadMaterial("materials/refractive_index/Au_44nm_Rosenblatt.csv")
Au_44nm_Rosenblatt(λ) = Au_44nm_Rosenblatt_estimator(λ)

Au_Babar_estimator = LoadMaterial("materials/refractive_index/Au_Babar.csv")
Au_Babar(λ) = Au_Babar_estimator(λ)

Au_McPeak_estimator = LoadMaterial("materials/refractive_index/Au_McPeak.csv")
Au_McPeak(λ) = Au_McPeak_estimator(λ)

Au_OlmonEvaporated_estimator = LoadMaterial("materials/refractive_index/Au_OlmonEvaporated.csv")
Au_OlmonEvaporated(λ) = Au_OlmonEvaporated_estimator(λ)

Au_OlmonSingleChrystaline_estimator = LoadMaterial("materials/refractive_index/Au_OlmonSingleChrystaline.csv")
Au_OlmonSingleChrystaline(λ) = Au_OlmonSingleChrystaline_estimator(λ)

Au_OlmonTemplateStripped_estimator = LoadMaterial("materials/refractive_index/Au_OlmonTemplateStripped.csv")
Au_OlmonTemplateStripped(λ) = Au_OlmonTemplateStripped_estimator(λ)


Au_Johnson_estimator = LoadMaterial("materials/refractive_index/Au_Johnson.csv")
Au_Johnson(λ) = Au_Johnson_estimator(λ)

Pd_Johnson_estimator = LoadMaterial("materials/refractive_index/Pd_Johnson.csv")
Pd_Johnson(λ) = Pd_Johnson_estimator(λ)


Au_Werner_estimator = LoadMaterial("materials/refractive_index/Au_Werner.csv")
Au_Werner(λ) = Au_Werner_estimator(λ)

Pd_Werner_estimator = LoadMaterial("materials/refractive_index/Pd_Werner.csv")
Pd_Werner(λ) = Pd_Werner_estimator(λ)


Ag_estimator = LoadMaterial("materials/refractive_index/Ag.csv")
Ag(λ) = Ag_estimator(λ)

LiF_estimator = LoadMaterial("materials/refractive_index/LiF.csv")
LiF(λ) = LiF_estimator(λ)

Pd_Palm_2018_estimator = LoadMaterial("materials/refractive_index/Pd_Palm.csv")
Pd_Palm_2018(λ) = Pd_Palm_2018_estimator(λ)


# F. Downes eq.1
# Note that small wavelengths give imaginary refractive index, this gives unphysical  results
SiO2_core_Sellmeier(λ) = √( 1
    + 0.6961663*(λ*1e6)^2/((λ*1e6)^2-0.0684043)
    + 0.4079426*(λ*1e6)^2/((λ*1e6)^2-0.1162414)
    + 0.8974794*(λ*1e6)^2/((λ*1e6)^2-9.896161) )

# F. Downes eq. 2
# Note that small wavelengths give imaginary refractive index, this gives unphysical  results
SiO2_thinfilm_Ciprian(λ) = √( 1 + 1.1336*(λ*1e6)^2/((λ*1e6)^2-9.261e-2) )


Au_unloaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/unloaded_n.csv", "materials/palm_PdAu-alloys/unloaded_k.csv", 1)
Au_unloaded(λ) = Au_unloaded_estimator(λ)

Pd014_unloaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/unloaded_n.csv", "materials/palm_PdAu-alloys/unloaded_k.csv", 2)
Pd014_unloaded(λ) = Pd014_unloaded_estimator(λ)

Pd034_unloaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/unloaded_n.csv", "materials/palm_PdAu-alloys/unloaded_k.csv", 3)
Pd034_unloaded(λ) = Pd034_unloaded_estimator(λ)

Pd042_unloaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/unloaded_n.csv", "materials/palm_PdAu-alloys/unloaded_k.csv", 4)
Pd042_unloaded(λ) = Pd042_unloaded_estimator(λ)

Pd052_unloaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/unloaded_n.csv", "materials/palm_PdAu-alloys/unloaded_k.csv", 5)
Pd052_unloaded(λ) = Pd052_unloaded_estimator(λ)

Pd073_unloaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/unloaded_n.csv", "materials/palm_PdAu-alloys/unloaded_k.csv", 6)
Pd073_unloaded(λ) = Pd073_unloaded_estimator(λ)

Pd_unloaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/unloaded_n.csv", "materials/palm_PdAu-alloys/unloaded_k.csv", 7)
Pd_unloaded(λ) = Pd_unloaded_estimator(λ)


Pd034_loaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/loaded_n.csv", "materials/palm_PdAu-alloys/loaded_k.csv", 3)
Pd034_loaded(λ) = Pd034_loaded_estimator(λ)

Pd042_loaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/loaded_n.csv", "materials/palm_PdAu-alloys/loaded_k.csv", 4)
Pd042_loaded(λ) = Pd042_loaded_estimator(λ)

Pd052_loaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/loaded_n.csv", "materials/palm_PdAu-alloys/loaded_k.csv", 5)
Pd052_loaded(λ) = Pd052_loaded_estimator(λ)

Pd073_loaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/loaded_n.csv", "materials/palm_PdAu-alloys/loaded_k.csv", 6)
Pd073_loaded(λ) = Pd073_loaded_estimator(λ)

Pd_loaded_estimator = LoadMaterial("materials/palm_PdAu-alloys/loaded_n.csv", "materials/palm_PdAu-alloys/loaded_k.csv", 7)
Pd_loaded(λ) = Pd_loaded_estimator(λ)


end # spesifics



end # materials