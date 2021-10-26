module utilities

export c_0, ϵ_0, peak_width

c_0 = 299792458
ϵ_0 = 8.854187812813e-12

function peak_width(data::Vector{<:Real}, step::T, inverted::Bool)::T where T<:Real
    data .-= minimum(data)
    if inverted
        data = maximum(data) .- data
    end

    target = maximum(data) / ℯ

    previous::Real = 0
    left::UInt = 0
    right::UInt = 0
    for (i, current) in zip(1:1:length(data), data)
        if left == 0
            if (previous < target) & (current >= target)
                left = i
            end
        elseif right == 0
            if (previous > target) & (current <= target)
                right = i
            end
        else
            break
        end
        previous = current
    end

    (right - left) * step
end

end # utilities
