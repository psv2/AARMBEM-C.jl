# vdW interactions in AARMBEM-C.jl

!!! warning

    vdW interaction free energy calculations have not been tested *at
    all* in AARMBEM-C.jl. The code for computing vdW interaction free
    energies should be the same as in AARMBEM.jl in the absence of
    internuclear couplings between separate nuclear bodies, but does
    change in the presence of such couplings; neither case has been
    tested in AARMBEM-C.jl. Therefore, users should use extreme
    caution with this code; no guarantees can be made about its
    stability or validity. This could of course change for the better
    with future releases.