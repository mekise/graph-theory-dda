function greensfun(a, b, ω, J::Stdd; imagshift=1E-23)
    k = k0(ω, J)
    dim = length(a)
    # r = norm(a-b) + imagshift*1im
    r = norm(a-b) + imagshift
    if dim == 1
        result = exp(1im*k*r)/(2im*k)
    elseif dim == 2
        result = (1/(4im))*hankelh1(0, k*r)
    elseif dim == 3
        result = exp(1im*k*r)/(4π*r)
    end
    return result
end