module ema

using Polynomials
using GeometricalPredicates

using ..analyticalmaterials

export n_MaxwellGarnett, n_Bruggeman


function n_MaxwellGarnett(ϵ_a::Number, ϵ_b::Number, f_A::Real)::Number
    # Hans Arwin, eq. 2.23
    # Coated spheres model
    # f_A is the fill factor

    ϵ_total = ϵ_b * (ϵ_a + 2*ϵ_b + 2*f_A*(ϵ_a - ϵ_b)) / (ϵ_a + 2*ϵ_b - f_A*(ϵ_a - ϵ_b))
    n(ϵ_total)
end

function n_MaxwellGarnett(ϵ_H::Number, ϵ_is::Vector{<:Number}, f_is::Vector{<:Real})::Number
    Σ_i = sum(f_i * (ϵ_i - ϵ_H) / (ϵ_i + 2*ϵ_H) for (ϵ_i, f_i) in zip(ϵ_is, f_is))

    ϵ_total = ϵ_H * (2*Σ_i + 1) / (1 - Σ_i)
    n(ϵ_total)
end


function n_Bruggeman(ϵ_a::Number, ϵ_b::Number, f_A::Real)::Number
    # Hans Arwin, eq. 2.26
    # ϵ determination from Jansson & Arwin

    f_B::Real = 1 - f_A
    u::Number = (3*f_A - 1)*ϵ_a + (3*f_B - 1)*ϵ_b
    ϵ_total_1 = (u + sqrt(u^2 + 8*ϵ_a*ϵ_b)) / 4
    ϵ_total_2 = (u - sqrt(u^2 + 8*ϵ_a*ϵ_b)) / 4


    if select_Wiener_circle(ϵ_total_1, ϵ_a, ϵ_b)
        return n(ϵ_total_1)
    elseif select_Wiener_circle(ϵ_total_2, ϵ_a, ϵ_b)
        return n(ϵ_total_2)
    else
        @error "No valid solution from Wiener limits, returning best match."
        return imag(ϵ_total_1) > imag(ϵ_total_2) ? ϵ_total_1 : ϵ_total_2
    end
end

function n_Bruggeman(ϵ_a::Number, ϵ_b::Number, ϵ_c::Number, f_a::Real, f_b::Real, f_c::Real)::Number
    # Solved by hand in journal page 37 from Hans Arwin

    k_0 = ϵ_a * ϵ_b * ϵ_c * (f_a + f_b + f_c)
    k_1 = 0
    k_2 = 0
    k_3 = 0

    for (f_i, ϵ_i, ϵ_2, ϵ_3) in ((f_a, ϵ_a, ϵ_b, ϵ_c), (f_b, ϵ_b, ϵ_a, ϵ_c), (f_c, ϵ_c, ϵ_b, ϵ_b))
        k_1 += f_i * ϵ_i * ϵ_2 * 2 + f_i * ϵ_i * ϵ_3 * 2 - f_i * ϵ_2 * ϵ_3
        k_2 += f_i * ϵ_i * 4 - f_i * ϵ_2 * 2 - f_i * ϵ_3 * 2
        k_3 -= 4 * f_i
    end

    ϵ_total_1, ϵ_total_2, ϵ_total_3 = roots(Polynomial([k_0, k_1, k_2, k_3]))

    for ϵ_potential in (ϵ_total_1, ϵ_total_2, ϵ_total_3)
        if select_Wiener_polygon(ϵ_potential, [ϵ_a, ϵ_b, ϵ_c])
            return ϵ_potential
        elseif (select_Wiener_circle(ϵ_potential, ϵ_a, ϵ_b) | select_Wiener_circle(ϵ_potential, ϵ_a, ϵ_c) | select_Wiener_circle(ϵ_potential, ϵ_b, ϵ_c))
            return ϵ_potential
        else
            @error "No valid solution from Wiener limits, returning best matach."
            return 0
        end
    end
end


function select_Wiener_circle(ϵ_s::Number, ϵ_1::Number, ϵ_2::Number)::Bool
    # Check if on line with origo, Jansson & Arwin algorithm

    if (real(ϵ_1) * imag(ϵ_2) - real(ϵ_2) * imag(ϵ_1)) < 1e-8 * abs(ϵ_1)
        function w(z::Number)::Number
            z * conj(ϵ_2 - ϵ_1) / abs(ϵ_2 - ϵ_1)
        end
        w_s, w_1, w_2 = w(ϵ_s), w(ϵ_1), w(ϵ_2)

        if real(w_2) < real(w_1)
            w_s, w_1, w_2 = w_s * -1, w_1 * -1, w_2 * -1
        end

        return (imag(w_s) + w_1 ≈ w_1) & (real(w_1) <= real(w_s) + 1e-8) & (real(w_s) <= real(w_2) + 1e-8)
    end

    z_0 = ϵ_1 * ϵ_2 * (conj(ϵ_1) - conj(ϵ_2)) / (conj(ϵ_1) * ϵ_2 - ϵ_1 * conj(ϵ_2))

    function ζ(z::Number)::Number
        (z - z_0) * conj(ϵ_2 - ϵ_1) / (abs(z_0) * abs(ϵ_2 - ϵ_1))
    end
    ζ_s, ζ_1, ζ_2 = ζ(ϵ_s), ζ(ϵ_1), ζ(ϵ_2)

    if imag((ζ_1 + ζ_2)/2) < 0
        ζ_s, ζ_1, ζ_2 = ζ_s * -1, ζ_1 * -1, ζ_2 * -1
    end

    return (abs(ζ_s) <= 1 + 1e-8) & (imag(ζ_s) >= imag((ζ_1 + ζ_2)/2) + 1e-8)
end

function select_Wiener_polygon(ϵ_s::Number, ϵ_v::Vector{<:Number})::Bool
    # Geometrical limits from Jansson & Arwin, solution with modern tools

    polygonpoints = Vector{Point2D}()
    for ϵ in ϵ_v
        push!(polygonpoints, Point(real(ϵ), imag(ϵ)))
    end
    poly = Polygon(polygonpoints...)

    inpolygon(poly, Point(real(ϵ_s), imag(ϵ_s)))
end

end # ema