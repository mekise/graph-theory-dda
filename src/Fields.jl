function incfield(ϕinput, scattpos, alphas, ω, J::Stdd; returnfield=false, imagshift=1E-23)
    M = intmatrix(scattpos, alphas, ω, J; imagshift=imagshift)
    ϕinc_tilde = M \ ϕinput
    if returnfield
        return ϕinc_tilde./alphas
    else
        return ϕinc_tilde
    end
end

function totfield(x, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    ϕinc_tilde = incfield(ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)
    n = size(scattpos, 1)
    ϕtot = 0. + 0im
    for i in 1:n
        ϕtot += ϕinc_tilde[i]*greensfun(x, scattpos[i, :], ω, J; imagshift=imagshift)
    end
    return ϕtot
end

function totfieldpolar(x, ϕinput, scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    r, θ, ϕ = x[1], x[2], x[3]
    return totfield([r*cos(θ)*sin(ϕ), r*sin(θ)*sin(ϕ), r*cos(ϕ)], ϕinput, scattpos, alphas, ω, J; imagshift=imagshift)
end