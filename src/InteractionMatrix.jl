function intmatrix(scattpos, alphas, ω, J::Stdd; imagshift=1E-23)
    n = size(scattpos, 1)
    M = zeros(ComplexF64, (n, n))
    for i in 1:n
        for j in 1:n
            M[i, j] = -greensfun(scattpos[i, :], scattpos[j, :], ω, J; imagshift=imagshift)
        end
    end
    M[diagind(M)] .= 1 ./ alphas
    return M
end