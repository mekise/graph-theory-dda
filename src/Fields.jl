function incfield(ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false)
    M = intmatrix(scattpos, alphas, ω, J; normalized=normalized)
    ϕinc = M \ ϕinput
    return ϕinc
end

function totfield(x, ϕinput, scattpos, alphas, ω, J::Stdd; normalized=false)
    ϕinc = incfield(ϕinput, scattpos, alphas, ω, J; normalized=normalized)
    n = size(scattpos, 1)
    dim = size(scattpos, 2)
    ϕtot = 0. + 0im
    for i in 1:n
        ϕtot += alphas[i]*ϕinc[i]*greensfun(x, scattpos[i, :], ω, J)
    end
    return ϕtot
end