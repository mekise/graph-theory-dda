abstract type Parameter end

struct Stdd{T<:Real} <: Parameter
    v::T
end

function k0(ω, J::Stdd)
    return ω / J.v
end

function α(ω, ω0, γ)
    return 1 / (ω0^2 - ω^2 - 1im*γ*ω)
end