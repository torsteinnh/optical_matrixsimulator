module targetfigures

using ..matrixcore
using ..fresnelltools

export scann_singleparameter, scan_plasmon_singleparameter, scan_plasmon_dualparameter, make_layered_tm_system, make_d3_system, make_d2_dualparameter_system


function scann_singleparameter(system::Function, scanrange::Vector{<:Number}, postprocessing::Function=v -> abs(v[2, 1])^2)::Tuple{Vector{<:Number}, Vector}
    # This function provides a unified interface for scanning structures for different design parameters.
    # It akes a parameter "system", a function returning (for instance) a total scattering matrix, a scans it for all elements in the scanrange and a postprocessing expression.
    # Returns the scanrange and the output data.
    # The default postprocessing stage takes a scattering matrix and returns the power reflection coefficient.

    output = Vector(undef, length(scanrange))

    for i in 1:length(scanrange)
        instance = system(scanrange[i])
        output[i] = postprocessing(instance)
    end

    scanrange, output
end


function plasmon_minima(system::Function, λs::Vector{Float64}, predipps::Int64=1)::Tuple{Float64, Float64}
    # This function finds the plasmon minima, as well as the relevant frequency at which it appears.
    # It takes as its arguments a function returning a scattering matrix as a function of a wavelength, a range of wavelengths and a number of predipps.
    # The predipps number is the number of local maxima the function can encounter before searching for a global minima.
    # This is useful for situations in which the plasmon dip is weaker than the high-frequency (low wavelength) absorption of the material.
    # Note however that, if set too high, the algorithm might measure a local minima pas the plasmon dip instead.
    # The function returns the plasmon minima power reflection and the wavelength at which it was found.

    minλ = λs[1]
    minr = abs(system(minλ)[2, 1])^2
    initialrise = predipps
    thisr = 1

    for λ in λs
        lastr = thisr
        thisr = abs(system(λ)[2, 1])^2

        if (initialrise >= 1) && (thisr < lastr)
            initialrise -= 1
            minr = thisr
            minλ = λ
        end

        if thisr < minr
            minr = thisr
            minλ = λ
        end
    end

    minr, minλ
end


function scan_minima(xs::Vector{<:Number}, ys::Vector{<:Number})::Tuple{Number, Number}
    # An optimized function for finding the minimum of a curve.

    @assert(length(xs) == length(ys))
    @assert(length(xs) > 0)

    @inbounds begin
        xmin = xs[1]
        ymin = ys[1]

        for i in 1:length(xs)
            if ys[i] < ymin
                xmin = xs[i]
                ymin = ys[i]
            end
        end
    end

    xmin, ymin
end


function scan_plasmon_singleparameter(system::Function, scanrange::Vector{<:Number}, λs::Vector{Float64}, predipps::Int64=1)::Tuple{Vector{<:Number}, Vector{Float64}, Vector{Float64}}
    # A tool for running a single parameter plasmon scan.
    # Takes the following arguments:
    #   system: A function of the scanparameter returning a function on the form λ -> S.
    #   scanrange: A range of scan parameters.
    #   λs: A range of wavelengths.
    #   predipps: The number of local maxima in the reflection spectrum to scan through before finding the plasmon minima.
    # Returns the scanrange, the minimum reflection coefficients and the relevant wavelengths.
    # Functions similarly to the single thickness scann one-off written in demos/thickness.jl.

    plasmonfinder(subsystem) = plasmon_minima(subsystem, λs, predipps)
    _, output = scann_singleparameter(system, scanrange, plasmonfinder)

    scanrange, [x[1] for x in output], [x[2] for x in output]
end


function scan_plasmon_dualparameter(system::Function, parameter1s::Vector{<:Number}, parameter2s::Vector{<:Number}, λs::Vector{Float64}, predipps::Int64=1)Tuple{Vector{<:Number}, Vector{<:Number}, Vector{<:Number}}
    # A tool for scanning for plasmon peaks across two parameters.
    # The function runs the single parameter scan on parameter 2 for each of parameter 1, it then plots the optimal parameter 2 as a function of parameter 1.
    
    minRP2s = Vector(undef, length(parameter1s))
    minλP2s = Vector(undef, length(parameter1s))
    
    for i in 1:length(parameter1s)
        p1 = parameter1s[i]

        _, reflections, wavelengths = scan_plasmon_singleparameter(system(p1), parameter2s, λs, predipps)
        minRP2, _ = scan_minima(parameter2s, reflections)
        minλP2, _ = scan_minima(parameter2s, wavelengths) # TODO this might be a suboptimal target function, try to find something better based on the plataus observed.

        minRP2s[i] = minRP2
        minλP2s[i] = minλP2
    end

    parameter1s, minRP2s, minλP2s
end



function make_layered_tm_system(ns::Vector{Function}, ds::Vector{Float64}, θi::Float64)::Function
    # A tool for generating λ -> S functions for arbitrary layered systems.

    @assert(length(ns) == length(ds))

    function S(λ)
        θw = θi

        elements = Vector(undef, 2 * length(ns) - 2)

        for i in 1:(length(ns) - 1)
            bulk = FresnellSlab(ns[i](λ), 2*π/λ, ds[i], θw)
            _, border, θw = FresnellBoundrary(ns[i](λ), ns[i + 1](λ), θw)

            elements[2*i - 1] = bulk
            elements[2*i] = border
        end

        CascadeScattering(elements)
    end

    S
end


function make_d3_system(n1::Function, n2::Function, n3::Function, n4::Function, n5::Function, d1::Float64, d2::Float64, d4::Float64, d5::Float64, θi::Float64)::Function
    function tm_singleparameter_system(d3::Float64)::Function
        make_layered_tm_system([n1, n2, n3, n4, n5], [d1, d2, d3, d4, d5], θi)
    end

    tm_singleparameter_system
end
function make_d2_dualparameter_system(n1::Function, n2::Function, n3::Function, n4::Function, n5::Function, d1::Float64, d4::Float64, d5::Float64, θi::Float64)::Function
    function tm_dualparameter_system(d2::Float64)::Function
        make_d3_system(n1, n2, n3, n4, n5, d1, d2, d4, d5, θi)
    end

    tm_dualparameter_system
end


end # targetfigures