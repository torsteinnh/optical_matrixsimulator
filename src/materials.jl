module materials

using Interpolations
using DelimitedFiles


export LoadMaterial, spesifics


function LoadMaterial(path)
    material = readdlm(path, ',', header=true)[1][:, 1:3]

    n_data = [material[:, 1] .* 1e-6, material[:, 2]]
    k_data = [material[:, 1] .* 1e-6, material[:, 3]]
    n_real = LinearInterpolation(n_data[1], n_data[2])
    n_complex = LinearInterpolation(k_data[1], k_data[2])

    n_total(λ) = n_real(λ) - abs(n_complex(λ)) * 1im
    n_total
end


module spesifics

import ..LoadMaterial
using ...analyticalmaterials


export Au, Pd, Ag, SiO2_core_Sellmeier, SiO2_thinfilm_Ciprian


Au_estimator = LoadMaterial("materials/Au.csv")
Au(λ) = Au_estimator(λ)


Pd_estimator = LoadMaterial("materials/Pd.csv")
Pd(λ) = Pd_estimator(λ)


Ag_estimator = LoadMaterial("materials/Ag.csv")
Ag(λ) = Ag_estimator(λ)


# F. Downes eq.1
# Note that small wavelengths give imaginary refractive index, this gives unphysical  results
SiO2_core_Sellmeier(λ) = √( 1
    + 0.6961663*(λ*1e6)^2/((λ*1e6)^2-0.0684043)
    + 0.4079426*(λ*1e6)^2/((λ*1e6)^2-0.1162414)
    + 0.8974794*(λ*1e6)^2/((λ*1e6)^2-9.896161) )

# F. Downes eq. 2
# Note that small wavelengths give imaginary refractive index, this gives unphysical  results
SiO2_thinfilm_Ciprian(λ) = √( 1 + 1.1336*(λ*1e6)^2/((λ*1e6)^2-9.261e-2) )


end # spesifics



end # materials