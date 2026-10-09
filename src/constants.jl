# XSTAR conventions and numerical constants used by more than one file. The physical constants
# (k, h, c, the Rydberg, ...) are in physical.jl. Numbers that
# belong to a single fit or rate are named in the file that uses them.

# Conversions to the precisions of ucalc's arithmetic. They narrow or widen the machine floats the
# database and the tests use and leave other numbers (dual numbers for automatic differentiation,
# BigFloat) alone, so that the rates can be differentiated.
to_single(x::Float64) = Float32(x)
to_single(x) = x
to_double(x::Union{Float16, Float32}) = Float64(x)
to_double(x) = x

const cx_unit = 1e-9                  # charge-exchange rates are fitted in 10⁻⁹ cm³ s⁻¹
const tiny = 1e-24                    # guard against division by zero and zero wavelengths

# exp() with the argument clamped, as in XSTAR's expo()
const expo_limit = 60.0
expo(x) = exp(clamp(x, -expo_limit, expo_limit))

# int() of Fortran for a default (32-bit) integer: a value out of range is not an error (the arm64 conversion saturates, as here)
fortran_int(x) = trunc(Int, clamp(x, Float64(typemin(Int32)), Float64(typemax(Int32))))

const Mb = 1e-18                      # cm² per Mb (cross sections are tabulated in Mb)
