# AARMBEM-C.jl

AARMBEM-C.jl is a code to compute electromagnetic (EM) interactions,
particularly fluctuational EM interactions, and linear phonon
propagation at the nanoscale, via the retarded many-body (RMB)
framework. In particular, the RMB framework computes fluctuational EM
interactions among material bodies that may scale from single atoms or
small compact molecules to infinitely long atomic-scale wires or
sheets, in vacuum or in the presence of macroscopic
bodies. AARMBEM-C.jl is a computational implementation of the RMB
framework, extended to account for the possibility of propagation of
phonons (which may transport momentum and energy) between material
bodies.

!!! note

    Effectively, AARMBEM-C.jl is a modification of the implementation
    of AARMBEM.jl to account for the possibility of phonon
    transmission between material bodies.

Currently, two codes are available. One computes van der Waals (vdW)
interaction energies in such systems, while the other computes heat
transfer powers in such systems. Thus, the focus is on fluctuational
phenomena. However, the theory and the API are general enough that
others may extend this code to perform computations of deterministic
EM phenomena, like absorbed or scattered powers from a specified
incident field/polarization source, local densities of states, et
cetera, with the inclusion of phonon transport.

!!! warning

    AARMBEM-C.jl has *not* been tested anywhere close to as
    extensively as AARMBEM.jl. Many of the functions may be broken or
    have unforeseen behaviors, so users should use a high degree of
    caution when using AARMBEM-C.jl and interpreting the results. This
    warning also applies to the documentation: there may be many
    mistakes and issues, and less care has been taken to prepare it.

## Name

The name "AARMBEM-C" is an acronym, expanded to "Ab-initio Atomistic
Retarded Many-Body Electromagnetics at the Mesoscale with
Conduction". The first part, "AARMBEM", may be pronounced as
"aarambham", identical to a word common in Indian languages meaning
"beginning" or "ab-initio"; the second part is simply the letter "C"..

## Installation

In the future, AARMBEM-C.jl may be properly packaged as a module that
may be installed and used. For now, the simplest way to install it is
to download and locally save the full directory of code from GitHub,
and running the code from other directories pointing to this one.