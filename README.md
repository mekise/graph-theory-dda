# Graph theory approach to exceptional points in wave scattering

This is the code I developed to solve wave scattering systems built on discrete dipole approximation. The code can handle N scatterer systems in arbitrary settings, but shines for cyclic polygon settings. These can be automatically generated with the sole number of scatterers and radius of the polygon inscription (see the paper for an accurate description of the system).

### Independently from the system analyzed, the package can (but is not limited to)
- generate and solve the non-linear system to find N-th order exceptional points
- show coalescence of eigenvalues and eigenvectors. The latter in the form of vanishing global Euclidean distance
- evaluate the total field over a sweep in the position parameter r
- evaluate the power output over a sweep in the position parameter r
- find the best imaginary shift in the resonant frequency to tune gain/loss of scatterers

### Authors
S. Scali, J. Anders, S. A. R. Horsley

### Abstract
In this paper, we apply graph-theory tools to a system of scatterers in the picture of discrete dipole approximation to better understand weak and strong coupling limits. In addition, in the presence of active scatterers, we show how we can generate arbitrary large order exceptional points (EP) in parameter space. This results in the sensitivity of the power output of the system being proportional to the inverse power of the probing parameter. To conclude, we interpret the parameter conditions for generating the N-th order EPs in a graph-theory picture.
