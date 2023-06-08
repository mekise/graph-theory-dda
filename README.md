# Graph theory approach to exceptional points in wave scattering
This is the code I developed to solve wave scattering systems built on the discrete dipole approximation. The code can handle N scatterer systems in arbitrary settings, but shines for cyclic polygon settings. These can be automatically generated with the sole number of scatterers and radius of the polygon inscription (see the paper for an accurate description of the system). These settings are particularly suitable for an efficient search of EPs.

### Independently from the system analyzed, the package can (but is not limited to)
- generate and solve the non-linear system to find N-th order exceptional points
- show coalescence of eigenvalues and eigenvectors. The latter in the form of vanishing global Euclidean distance
- evaluate the total field over a sweep in the position parameter r
- evaluate the power output over a sweep in the position parameter r
- perform a sweep in imaginary shift of the design frequency to find the preferred dissipation balance in the set of scatterers

### Paper
https://arxiv.org/abs/2301.08257

### Authors
S. Scali, J. Anders, S. A. R. Horsley

### Abstract
In this paper, we use graph theory to solve wave scattering problems in the discrete dipole approximation. As a key result of this work, in the presence of active scatterers, we present a systematic method to find arbitrary large-order zero eigenvalue exceptional points (EPs). This is achieved by solving a set of non-linear equations that we interpret, in a graph theory picture, as vanishing sums of scattering events. We then show how the total field of the system responds to parameter perturbations at the EP. Finally, we investigate the sensitivity of the power output to imaginary perturbation in the design frequency. This perturbation can be employed to trade sensitivity for a different dissipation balance of the system. The purpose of the results of this paper is manifold. On the one hand, we aim to shed light on the link between graph theory and wave scattering. On the other hand, the results of this paper find application in all those settings where zero eigenvalue EPs play a unique role like in coherent perfect absorption (CPA) structures.
