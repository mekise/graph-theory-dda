function intmatrix(scattpos, alphas, ω, J::Stdd; normalized=false, imagshift=1E-12)
    n = size(scattpos, 1)
    M = zeros(ComplexF64, (n, n))
    for i in 1:n
        for j in 1:n
            M[i, j] = -alphas[i]*greensfun(scattpos[i, :], scattpos[j, :], ω, J; imagshift=imagshift)
        end
    end
    M[diagind(M)] .= 1.0
    if !normalized
        for i in 1:n
            M[i, :] /= alphas[i]
        end
    end
    return M
end