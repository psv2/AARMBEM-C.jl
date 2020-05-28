# Matrix storage conventions

This page goes over conventions for matrix storage in
AARMBEM-C.jl. Additionally, there are subtle differences in how matrices
are stored when computing vdW interaction free energy integrands
versus when computing thermal radiation spectra, so these are
discussed too.

## Matrix element ordering

The matrices ``K^{(0)}_{\mathrm{I}}`` and ``G^{(0)}`` are stored separately
for each individual molecular body labeled ``n``. Such a matrix ``A``
for an individual molecular body with ``N`` atoms, with elements
``A_{pi,qj}`` for atoms ``p, q \in \{1, 2, \ldots, N\}`` and Cartesian
direction indices ``i, j \in \{x, y, z\}`` is stored in the form
```math
A = \begin{bmatrix}
A_{1x,1x} & A_{1x,1y} & A_{1x,1z} & A_{1x,2x} & A_{1x,2y} & A_{1x,2z} & \ldots & A_{1x,Nx} & A_{1x,Ny} & A_{1x,Nz} \\
A_{1y,1x} & A_{1y,1y} & A_{1y,1z} & A_{1y,2x} & A_{1y,2y} & A_{1y,2z} & \ldots & A_{1y,Nx} & A_{1y,Ny} & A_{1y,Nz} \\
A_{1z,1x} & A_{1z,1y} & A_{1z,1z} & A_{1z,2x} & A_{1z,2y} & A_{1z,2z} & \ldots & A_{1z,Nx} & A_{1z,Ny} & A_{1z,Nz} \\
A_{2x,1x} & A_{2x,1y} & A_{2x,1z} & A_{2x,2x} & A_{2x,2y} & A_{2x,2z} & \ldots & A_{2x,Nx} & A_{2x,Ny} & A_{2x,Nz} \\
A_{2y,1x} & A_{2y,1y} & A_{2y,1z} & A_{2y,2x} & A_{2y,2y} & A_{2y,2z} & \ldots & A_{2y,Nx} & A_{2y,Ny} & A_{2y,Nz} \\
A_{2z,1x} & A_{2z,1y} & A_{2z,1z} & A_{2z,2x} & A_{2z,2y} & A_{2z,2z} & \ldots & A_{2z,Nx} & A_{2z,Ny} & A_{2z,Nz} \\
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
A_{Nx,1x} & A_{Nx,1y} & A_{Nx,1z} & A_{Nx,2x} & A_{Nx,2y} & A_{Nx,2z} & \ldots & A_{Nx,Nx} & A_{Nx,Ny} & A_{Nx,Nz} \\
A_{Ny,1x} & A_{Ny,1y} & A_{Ny,1z} & A_{Ny,2x} & A_{Ny,2y} & A_{Ny,2z} & \ldots & A_{Ny,Nx} & A_{Ny,Ny} & A_{Ny,Nz} \\
A_{Nz,1x} & A_{Nz,1y} & A_{Nz,1z} & A_{Nz,2x} & A_{Nz,2y} & A_{Nz,2z} & \ldots & A_{Nz,Nx} & A_{Nz,Ny} & A_{Nz,Nz}
\end{bmatrix}
```
in AARMBEM-C.jl, so this is a ``3N \times 3N`` matrix.

The intermolecular internuclear coupling and Green's function matrices
``\Delta K_{\mathrm{I}}`` and ``G^{\mathrm{mac}}`` must be stored for
all molecular bodies together. Such a matrix ``X`` for the whole
system has matrix blocks ``X_{mn}`` for molecules ``m, n \in \{1, 2,
\ldots, N_{\mathrm{mol}}\}``. The actual storage in AARMBEM-C.jl is as
a full matrix, but representing it compactly in this page is
difficult. Instead, it will be represented in this page in block form
as
```math
X = \begin{bmatrix}
X_{1,1} & X_{1,2} & \ldots & X_{1,N_{\mathrm{mol}}} \\
X_{2,1} & X_{2,2} & \ldots & X_{2,N_{\mathrm{mol}}} \\
\vdots & \vdots & \ddots & \vdots \\
X_{N_{\mathrm{mol}},1} & X_{N_{\mathrm{mol}},2} & \ldots & X_{N_{\mathrm{mol}},N_{\mathrm{mol}}}
\end{bmatrix}
```
where each block is expanded as
```math
X_{mn} = \begin{bmatrix}
(X_{mn})_{1x,1x} & (X_{mn})_{1x,1y} & (X_{mn})_{1x,1z} & (X_{mn})_{1x,2x} & (X_{mn})_{1x,2y} & (X_{mn})_{1x,2z} & \ldots & (X_{mn})_{1x,N_{n}x} & (X_{mn})_{1x,N_{n}y} & (X_{mn})_{1x,N_{n}z} \\
(X_{mn})_{1y,1x} & (X_{mn})_{1y,1y} & (X_{mn})_{1y,1z} & (X_{mn})_{1y,2x} & (X_{mn})_{1y,2y} & (X_{mn})_{1y,2z} & \ldots & (X_{mn})_{1y,N_{n}x} & (X_{mn})_{1y,N_{n}y} & (X_{mn})_{1y,N_{n}z} \\
(X_{mn})_{1z,1x} & (X_{mn})_{1z,1y} & (X_{mn})_{1z,1z} & (X_{mn})_{1z,2x} & (X_{mn})_{1z,2y} & (X_{mn})_{1z,2z} & \ldots & (X_{mn})_{1z,N_{n}x} & (X_{mn})_{1z,N_{n}y} & (X_{mn})_{1z,N_{n}z} \\
(X_{mn})_{2x,1x} & (X_{mn})_{2x,1y} & (X_{mn})_{2x,1z} & (X_{mn})_{2x,2x} & (X_{mn})_{2x,2y} & (X_{mn})_{2x,2z} & \ldots & (X_{mn})_{2x,N_{n}x} & (X_{mn})_{2x,N_{n}y} & (X_{mn})_{2x,N_{n}z} \\
(X_{mn})_{2y,1x} & (X_{mn})_{2y,1y} & (X_{mn})_{2y,1z} & (X_{mn})_{2y,2x} & (X_{mn})_{2y,2y} & (X_{mn})_{2y,2z} & \ldots & (X_{mn})_{2y,N_{n}x} & (X_{mn})_{2y,N_{n}y} & (X_{mn})_{2y,N_{n}z} \\
(X_{mn})_{2z,1x} & (X_{mn})_{2z,1y} & (X_{mn})_{2z,1z} & (X_{mn})_{2z,2x} & (X_{mn})_{2z,2y} & (X_{mn})_{2z,2z} & \ldots & (X_{mn})_{2z,N_{n}x} & (X_{mn})_{2z,N_{n}y} & (X_{mn})_{2z,N_{n}z} \\
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
(X_{mn})_{N_{m}x,1x} & (X_{mn})_{N_{m}x,1y} & (X_{mn})_{N_{m}x,1z} & (X_{mn})_{N_{m}x,2x} & (X_{mn})_{N_{m}x,2y} & (X_{mn})_{N_{m}x,2z} & \ldots & (X_{mn})_{N_{m}x,N_{n}x} & (X_{mn})_{N_{m}x,N_{n}y} & (X_{mn})_{N_{m}x,N_{n}z} \\
(X_{mn})_{N_{m}y,1x} & (X_{mn})_{N_{m}y,1y} & (X_{mn})_{N_{m}y,1z} & (X_{mn})_{N_{m}y,2x} & (X_{mn})_{N_{m}y,2y} & (X_{mn})_{N_{m}y,2z} & \ldots & (X_{mn})_{N_{m}y,N_{n}x} & (X_{mn})_{N_{m}y,N_{n}y} & (X_{mn})_{N_{m}y,N_{n}z} \\
(X_{mn})_{N_{m}z,1x} & (X_{mn})_{N_{m}z,1y} & (X_{mn})_{N_{m}z,1z} & (X_{mn})_{N_{m}z,2x} & (X_{mn})_{N_{m}z,2y} & (X_{mn})_{N_{m}z,2z} & \ldots & (X_{mn})_{N_{m}z,N_{n}x} & (X_{mn})_{N_{m}z,N_{n}y} & (X_{mn})_{N_{m}z,N_{n}z}
\end{bmatrix}
```
so these blocks are each ``3N_{m} \times 3N_{n}`` matrices where
molecular body ``m`` has ``N_{m}`` atoms and molecular body ``n`` has
``N_{n}`` atoms.

## Matrices in thermal energy transport powers

For thermal energy transport powers and vdW interaction free energies,
the total response matrix of the system is paramount. However, because
both electronic and nuclear oscillators contribute to the sources and
sinks of energy transport (as opposed to AARMBEM.jl where effectively
only electronic oscillators were sinks of energy transport due to the
lack of consideration of phonon conduction), the only full matrices
that make sense to store for each individual molecular body are
``G^{(0)}_{nn}`` and ``K^{(0)}_{\mathrm{I}n}``. The matrices
``G^{\mathrm{mac}}`` and ``\Delta K_{\mathrm{I}}`` are computed for
the full system, and from this the terms needed to compute the energy
transfer spectrum are constructed essentially on the fly, because so
many different matrices are required, and it is much more efficient to
perform multiplications of dense matrices with diagonal matrices
as-needed instead of densely storing such matrices at effectively
double the required size (in terms of memory). The computed quantities
can be derived by analytically computing the block inversion ``Y =
(Z^{(0)} + \Delta Z)^{-1}``, but for brevity, this will not be done in
this documentation.

Note that for thermal energy transport, ``G^{\mathrm{mac}}`` includes
the diagonal vacuum self-interaction elements ``\langle \vec{f}_{pi},
\mathbb{G}^{(0)} \vec{f}_{pi} \rangle`` only for evaluating the term
``\operatorname{asym}(P_{n} \Delta Z)`` in order to capture far-field
EM thermal emission contributions in vacuum, *not* for evaluating the
total response matrix ``Y``.

## Matrices in vdW interaction free energies

There are a few differences in the matrices used for computing vdW
interactions compared to those for computing thermal radiation. These
changes are for the sake of greater efficiency and numerical
accuracy. The expressions below are quoted for generic ``\omega`` but
evaluated specifically for ``\omega = \mathrm{i}\xi``. In particular,
the diagonal matrices ``Y^{(0)}_{\mathrm{e}n} = (K_{\mathrm{e}n} -
\mathrm{i}\omega B_{\mathrm{e}n} - \omega^{2} M_{\mathrm{e}n})^{-1}``
and ``Y^{(0)}_{\mathrm{I}n} = (K_{\mathrm{e}n} - \mathrm{i}\omega
B_{\mathrm{I}n} - \omega^{2} M_{\mathrm{I}n})^{-1}`` are frequently
used as preconditioners for more efficient and stable
computation. Additionally, as ``G^{\mathrm{mac}}`` is only used to
construct ``Y`` and ``Y_{\infty}``, it never includes the diagonal
vacuum self-interaction elements ``\langle \vec{f}_{pi},
\mathbb{G}^{(0)} \vec{f}_{pi} \rangle``.