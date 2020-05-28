# Morse potential for intermolecular internuclear couplings

AARMBEM-C.jl allows for the possibility of coherent thermal phonon
transport among molecular bodies when computing vdW interactions and
thermal energy transport, and these couplings are computed from the
Morse potential. In particular, a pair of atoms has equilibrium
separation ``a_{0}``, dissociation energy ``E_{\infty}``, and
equilibrium spring constant ``k_{0}``, then the spring constant at
separation ``r`` is
```math
k(r) = \frac{\sqrt{2k_{0} E_{\infty}}}{r - a_{0}} \left(\exp\left(-\sqrt{\frac{2k_{0}}{E_{\infty}}} (r - a_{0})\right) - \exp\left(-\sqrt{\frac{k_{0}}{2E_{\infty}}} (r - a_{0})\right)\right)
```
which is always negative and is equal to ``-k_{0}`` for ``r = a_{0}``.

The units are consistent with SI. In particular, ``a_{0}`` should be
given in meters, ``k_{0}`` in newtons per meter, ``E_{\infty}`` in
joules, and ``r`` in meters.