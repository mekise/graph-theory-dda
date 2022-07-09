include("./ScatterersDDA.jl")
using .ScatterersDDA

ω = 1.
J = Stdd(2);
scattpos = [[1. 2. 3.]
            [2. 3. 4.]
            [3. 3. 3.]
            [2. 1. 3.]];
alphas = [(1. + 2im)
          (1. - 2im)
          (2. + 1im)
          (1. - 3im)];
ϕinput = [(1. +1im), 2., 3., 2.];

incfield(ϕinput, scattpos, alphas, ω, J)
