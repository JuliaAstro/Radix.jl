# Temporary arrays of the rates that are not allocated at each call: the conversion of a table to `Float64` on access, and the work arrays of the integrals over the energy grid that each task of Julia keeps.

"""
    as_type(T, v)

The array `v` with its elements converted to `T` on access, with no copy: `T.(v)` for the rates whose tables are stored as `Float32` (the records of the database) and only read. `v` itself if it already has the type `T`.
`as_float64(v)` is `as_type(Float64, v)`.
"""
struct Converted{T, S, N, V<:AbstractArray{S, N}} <: AbstractArray{T, N}
    parent::V
end
Base.size(v::Converted) = size(v.parent)
Base.IndexStyle(::Type{<:Converted}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::Converted{T}, i::Int) where {T} = T(v.parent[i])
as_type(::Type{T}, v::AbstractArray{T}) where {T} = v
as_type(::Type{T}, v::AbstractArray) where {T} = Converted{T, eltype(v), ndims(v), typeof(v)}(v)
as_float64(v) = as_type(Float64, v)

"""
    scratch(key, n, T=Float64)

A work array of `n` elements of the type `T` that the current task of Julia keeps under the name `key` (its memory grows as needed, and its contents are those of the last use: not zeroed). Two arrays in use at the same time need different keys.
"""
function scratch(key::Symbol, n::Integer, ::Type{T}=Float64) where {T}
    store = get!(Dict{Symbol, Any}, task_local_storage(), :RadixScratch)::Dict{Symbol, Any}
    v = get(store, key, nothing)
    if !(v isa Vector{T})
        v = Vector{T}(undef, n)
        store[key] = v
    end
    length(v) == n || resize!(v, n)
    v
end

"""
    scaled_table(xs, f, key; floor=nothing, cap=nothing)

`max.(min.(xs*f, cap), floor)` (the limits that are given) of the table `xs`, in the work array `key` of the task if the numbers are `Float64` (valid until the next use of `key`), a new array for other numbers (dual numbers).
"""
function scaled_table(xs, f, key::Symbol; floor=nothing, cap=nothing)
    T = promote_type(eltype(xs), typeof(f))
    out = T === Float64 ? scratch(key, length(xs)) : Vector{T}(undef, length(xs))
    for i in eachindex(xs)
        v = xs[i]*f
        cap === nothing || (v = min(v, cap))
        floor === nothing || (v = max(v, floor))
        out[i] = v
    end
    out
end
