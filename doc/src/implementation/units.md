# Units in AARMBEM-C.jl

AARMBEM-C.jl uses a "Lorentz-Heaviside SI" unit system. This means that
the mass, length, and time units are the standard SI units of
kilograms, meters, and seconds. However, EM units are defined by the
Lorentz-Heaviside convention for Maxwell's equations, namely
```math
\nabla \cdot \vec{E} = \rho \\
\nabla \cdot \vec{B} = 0 \\
\nabla \times \vec{E} = -\frac{1}{c} \frac{\partial \vec{B}}{\partial t} \\
\nabla \times \vec{B} = \frac{1}{c} \left(\vec{J} + \frac{\partial \vec{E}}{\partial t}\right)
```
such that ``\vec{E}`` and ``\vec{B}`` have the same units as each
other. In material media, they also have the same units as ``\vec{D}``
and ``\vec{H}``. This is trivially true for the magnetic field, with
``\vec{B} = \vec{H}``, as AARMBEM-C.jl currently does not support
consideration of nontrivial magnetic susceptibilities. For the
electric field, this means ``\vec{D} = \vec{E} + \vec{P}``, where in
the frequency domain ``P_{i} (\vec{x}) = \sum_{j} \int V_{ij}
(\vec{x}, \vec{x}') E_{j} (\vec{x}')~\mathrm{d}^{3} x'``, so the
frequency domain kernel ``V_{ij} (\vec{x}, \vec{x}')`` has units of
inverse spatial volume. For a local homogeneous isotropic medium with
``V_{ij} (\vec{x}, \vec{x}') = \chi \Theta(\vec{x} \in V) \delta_{ij}
\delta^{3} (\vec{x} - \vec{x}')``, the susceptibility ``\chi`` is
dimensionless and is related to the permittivity ``\epsilon = 1 +
\chi``: this aligns with the SI convention (without the factor of
``\epsilon_{0}``), *not* the Gaussian CGS convention. Additionally, in
the frequency domain, the vacuum Maxwell Green's function is defined
to satisfy
```math
\sum_{k} (\partial_{i} \partial_{k} G^{(0)}_{kj} (\vec{x}, \vec{x}') -
\partial_{k} \partial_{k} G^{(0)}_{ij} (\vec{x}, \vec{x}')) -
(\omega/c)^{2} G^{(0)}_{ij} (\vec{x}, \vec{x}') = (\omega/c)^{2}
\delta_{ij} \delta^{3} (\vec{x} - \vec{x}')
```
with the explicit outgoing solution ``G^{(0)}_{ij} (\vec{x}, \vec{x}')
= [\partial_{i} \partial_{j} + (\omega/c)^{2}
\delta_{ij}](e^{\mathrm{i}\omega\vert\vec{x} -
\vec{x}'\vert/c}/(4\pi\vert\vec{x} - \vec{x}'\vert))``, where
``\partial_{i} = \frac{\partial}{\partial x_{i}}``; this means
``G^{(0)}_{ij} (\vec{x}, \vec{x}')``, along with
``G^{\mathrm{mac}}_{ij}(\vec{x}, \vec{x}')`` for any macroscopic body,
have units in position space of inverse spatial volume.

## Units of input parameters

Various functions require evaluation at frequency ``\omega`` (whether
real or imaginary), positions ``\vec{r}_{p}`` and ``\vec{r}_{q}``, and
optionally the Bloch wavevector ``\vec{k}``. The frequencies
``\omega`` are consistently taken to be *angular* frequencies in units
of radians per second: the Fourier expansion in time in terms of
frequency is defined to be ``g(t) = \int_{-\infty}^{\infty}
e^{-\mathrm{i}\omega t}
\tilde{g}(\omega)~\frac{\mathrm{d}\omega}{2\pi}``, and the Fourier
transform from time to frequency is ``\tilde{g}(\omega) =
\int_{-\infty}^{\infty} e^{\mathrm{i}\omega t} g(t)~\mathrm{d}t``. The
components of wavevectors ``\vec{k}`` are consistently taken to be
*angular* wavevectors in units of radians per meter: the Fourier
expansion in position in terms of wavevector is defined to be
``f(\vec{x}) = \int e^{\mathrm{i}\vec{k} \cdot \vec{x}}
\tilde{f}(\vec{k})~\frac{\mathrm{d}^{d} k}{(2\pi)^{d}}``, and the
Fourier transform from position to wavevector is ``\tilde{f}(\vec{k})
= \int e^{-\mathrm{i}\vec{k} \cdot \vec{x}} f(\vec{x})~\mathrm{d}^{d}
x``. The positions simply use the standard SI unit of meters for all
vector components.

Beyond these, there are several key inputs required by
AARMBEM-C.jl. The masses ``M_{\mathrm{e}}`` and ``M_{\mathrm{I}}`` are
given simply in kilograms. The damping coefficients ``B_{\mathrm{e}}``
and ``B_{\mathrm{I}}`` are given simply in kilograms per second. The
spring constants ``K_{\mathrm{e}}``, ``K^{(0)}_{\mathrm{I}}``, and
``\Delta K_{\mathrm{I}}`` are given simply in kilograms per second
squared (or, equivalently, newtons per meter). Only the charges
``Q_{\mathrm{e}}`` must be treated with a little more care: they must
be given in Lorentz-Heaviside SI units, where the conversion formula
from standard SI may be written as ``q_{\mathrm{LHSI}} =
q_{\mathrm{SI}} / \sqrt{\epsilon_{0(\mathrm{SI})}}``, and the result
is specified in ``\mathrm{kg}^{1/2} \cdot \mathrm{m}^{3/2} \cdot
\mathrm{s}^{-1}``.

## Units of output polarizabilities and Green's functions

AARMBEM-C.jl frequently makes use of atomic polarizabilities as well as
Green's function components, whether as effective per-atom quantities
or as full matrices. In either case, polarizabilities are in units of
cubic meters (so multiplying by ``\epsilon_{0}``, with no additional
factors of ``4\pi``, recovers the SI values), and Green's functions
are in units of inverse cubic meters.