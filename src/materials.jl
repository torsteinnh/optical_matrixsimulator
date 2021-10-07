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
using ...drudemetals


export Au, Au_drude


Au_estimator = LoadMaterial("materials/Au.csv")
Au(λ) = Au_estimator(λ)

Au_drude(λ) = n_drude(44.2e6, 27.3e-15, λ)

end # spesifics



end # materials