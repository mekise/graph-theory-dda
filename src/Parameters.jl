abstract type Parameters end

struct Complete{T<:Real} <: Parameters
    ω::T
    v::T
end

function k0(J::Complete)
    return J.ω / J.v
end

function α(ω0, γ, J::Complete)
    return 1 / (ω0^2 - J.ω^2 - 1im*γ*J.ω)
end