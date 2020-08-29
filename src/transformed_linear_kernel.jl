
struct Transform{F <: Function}
    f::F
end

(t::Transform)(x::AV{<:Real}) = t.f.(x)
(t::Transform)(x::ColVecs) = ColVecs(hcat(t.f.(x)...)) # In slow method error :<
# (t::Transform)(x::ColVecs) = hcat(t.f.(x)...) # In slow method
# (t::Transform)(x::ColVecs) = ColVecs(mapslices(t.f, x.X, dims=1))

#=
mapslices implementation make no sense since it finally behaves like "slow" method if
    it follows "AD style" to get gradient.

So we just disable ColVecs's error
=#

Zygote.@adjoint function ColVecs(X::AbstractMatrix)
    back(Δ::NamedTuple) = (Δ.X,)
    back(Δ::AbstractMatrix) = (Δ,)
    function back(Δ::AbstractVector{<:AbstractVector{<:Real}})
        (hcat(Δ...),)
    end
    return ColVecs(X), back
end


"""
    LinearWithTransform{T<:Real} <: BaseKernel
The standardised linear kernel / dot-product kernel with an extra transform.
"""
struct LinearWithTransform{T <: Transform} <: Stheno.BaseKernel
    transform::T
end

# Binary methods
function ew(k::LinearWithTransform, x::AV{<:Real}, x′::AV{<:Real})
    x = k.transform(x)
    x′ = k.transform(x′)
    x .* x′
end
function pw(k::LinearWithTransform, x::AV{<:Real}, x′::AV{<:Real})
    x = k.transform(x)
    x′ = k.transform(x′)
    x .* x′'
end
function ew(k::LinearWithTransform, x::ColVecs, x′::ColVecs)
    x = k.transform(x)
    x′ = k.transform(x′)
    reshape(sum(x.X .* x′.X; dims=1), :)
end
function pw(k::LinearWithTransform, x::ColVecs, x′::ColVecs)
    x = k.transform(x)
    x′ = k.transform(x′)
    x.X' * x′.X
end

# Unary methods
function ew(k::LinearWithTransform, x::AV{<:Real})
    x = k.transform(x)
    x.^2
end
function pw(k::LinearWithTransform, x::AV{<:Real})
    x = k.transform(x)
    x .* x'
end
function ew(k::LinearWithTransform, x::ColVecs)
    x = k.transform(x)
    reshape(sum(abs2.(x.X); dims=1), :)
end
function pw(k::LinearWithTransform, x::ColVecs)
    x = k.transform(x)
    x.X' * x.X
end