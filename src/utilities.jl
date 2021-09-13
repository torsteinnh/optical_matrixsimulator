module utilities

export c_0

function Base.show(io::IO, num::Complex)
    # A simple overloading of the default complex number printing for polar coordinates
    show(io, abs(num))
    write(io, "ℯ^{")
    show(io, angle(num) / π)
    write(io, "π}")
end

c_0 = 299792458

end # utilities
