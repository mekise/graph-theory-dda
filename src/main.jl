include("./ScatterersDDA.jl")
using .ScatterersDDA

ω = 1.
J = Standard(2.);

a = [1., 2., 3.]
b = [2., 3., 2.]

G = G(a, b, ω, J)

print(G)