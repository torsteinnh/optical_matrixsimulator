module materials

using Interpolations
using DelimitedFiles

export LoadMaterial


function LoadMaterial(path)
    material = readdlm(path, ',', header=true)[1][:, 1:3]

    n_data = [material[:, 1] .* 1e-6, material[:, 2]]
    k_data = [material[:, 1] .* 1e-6, material[:, 3]]
    n_real = LinearInterpolation(n_data[1], n_data[2])
    n_complex = LinearInterpolation(k_data[1], k_data[2])

    n_total(λ) = n_real(λ) + n_complex(λ) * 1im
    n_total
end



end # materials