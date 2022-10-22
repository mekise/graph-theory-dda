module ScatterersDDA

using LinearAlgebra
using SpecialFunctions
using Random
using ForwardDiff
using Cubature
using Combinatorics

include("Misc.jl");
include("GreensFunction.jl");
include("InteractionMatrix.jl");
include("Fields.jl");
include("ComplexGradient.jl")
include("Power.jl")
include("Approximation.jl")

export Parameter, Stdd
export k0, α, greensfun, intmatrix, incfield, totfield, totfieldpolar, complexgrad, complexder
export integrand, powerout, poweroutexplicit, evalsumm, evalsummcorrected
export determinant, det_weak, permsingleinversion, weak4by4, weakincfield, weaktotfield, strong4by4, strongincfield, strongtotfield

end