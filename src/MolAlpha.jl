"""

    constructMolAlpha!(myMolData::OneMol_WithPhonons{FT},
                       myPeriodicData::PeriodicData{FT},
                       freq::Union{FT, Complex{FT}},
                       k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the vector of atomic polarizabilities `myMolData.alpha0` based on the
susceptibility of the molecule in isolation. `myMolData.Ye0ee` will incorporate the
electronic contributions to the response, and `myMolData.Ye0II` will incorporate most
aspects of the nuclear contributions to the response, mainly of use as preconditioners when
computing vdW interactions.

"""
function constructMolAlpha!(myMolData::OneMol_WithPhonons{FT},
                            myPeriodicData::PeriodicData{FT},
                            freq::Union{FT, Complex{FT}},
                            k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    alphainv = zeros(Complex{FT}, 3*myMolData.numatoms, 3*myMolData.numatoms);
    Ke3N = vecNto3N(myMolData.Ke);
    Qeinv3N = vecNto3N(1 ./ myMolData.Qe);

    KIk!(myMolData, myPeriodicData, k);
  
    setindex!(alphainv, myMolData.KIBloch,
              1:3*myMolData.numatoms, 1:3*myMolData.numatoms);
    addtodiag!(alphainv,
               Ke3N .- vecNto3N(1im.*freq.*myMolData.BI .+ freq^2 .* myMolData.MI));
    setindex!(alphainv, inv(alphainv),
              1:3*myMolData.numatoms, 1:3*myMolData.numatoms);
    mulMatDiagright!(alphainv, Ke3N);
    mulMatDiagleft!(alphainv, Ke3N);
    subtractfromdiag!(alphainv,
                      Ke3N .- vecNto3N(1im.*freq.*myMolData.Be .+ freq^2 .* myMolData.Me));
    mulMatDiagright!(alphainv, Qeinv3N);
    mulMatDiagleft!(alphainv, Qeinv3N);
    alpha = inv(alphainv);
    alphacarttr = Array{eltype(alpha), 2}(undef, myMolData.numatoms, myMolData.numatoms);
    @inbounds for qq=1:myMolData.numatoms
        @inbounds for pp=1:myMolData.numatoms
            alphacarttr[pp, qq] = tr(alpha[3*pp-2:3*pp, 3*qq-2:3*qq]);
        end
    end
  
    setindex!(myMolData.alpha0, abs.(vec(sum(alphacarttr, dims=2))), 1:myMolData.numatoms);
    lmul!(1/3, myMolData.alpha0);
    setindex!(myMolData.Ye0ee,
              1 ./ (myMolData.Ke .-
                    (1im .* freq .* myMolData.Be .+ freq^2 .* myMolData.Me)),
              1:myMolData.numatoms);
    setindex!(myMolData.Ye0II,
              1 ./ (myMolData.Ke .-
                    (1im .* freq .* myMolData.BI .+ freq^2 .* myMolData.MI)),
              1:myMolData.numatoms);
  
end


## TO DELETE: constructMolAlpha!(myMolData::OneMol_NoPhonons{FT}, ...)
## because it no longer makes sense for phonon conduction (if phonons
## are to be ignored, use the usual RMB code)
"""

    constructMolAlpha!(myMolData::OneMol_NoPhonons{FT},
                       myPeriodicData::PeriodicData{FT},
                       freq::Union{FT, Complex{FT}},
                       k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the inverse molecular susceptibility `myMolData.alphainv` in-place at frequency
`freq` and wavevector `k`, the latter of which must be a 3-element real vector. Without
phonons, `myMolData.alphainv` will be a vector, and the vectors `myMolData.alpha0` and
`myMolData.alphae` will be simply related to `myMolData.alphainv`.

"""
function constructMolAlpha!(myMolData::OneMol_NoPhonons{FT},
                            myPeriodicData::PeriodicData{FT},
                            freq::Union{FT, Complex{FT}},
                            k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    setindex!(myMolData.alphainv,
              (myMolData.Ke .- (1im.*freq.*myMolData.Be .+
                                freq^2 .* myMolData.Me))./(myMolData.Qe.^2),
              1:myMolData.numatoms);
    setindex!(myMolData.alphae, 1 ./ myMolData.alphainv, 1:myMolData.numatoms);
    setindex!(myMolData.alpha0, abs.(myMolData.alphae), 1:myMolData.numatoms);
    ## constructMolAlphaNoPhonons!(myMolData, myPeriodicData, freq, k)
  
end

"""

    KIk(myMolData::OneMol_WithPhonons{FT}, myPeriodicData::PeriodicData{FT},
        k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the `k`-dependence of the internuclear spring constant matrix for the molecule in
isolation, and return the result. Return `myMolData.KI` if the geometry has no periodic
dimensions.

See also: [`constructMolAlpha!`](@ref).

"""
function KIk(myMolData::OneMol_WithPhonons{FT}, myPeriodicData::PeriodicData{FT},
             k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    if (myPeriodicData.numdims < 1)
        return complex(myMolData.KI);
    end
    if (length(k) != 3)
        error("KIk(): k must have length 3");
    end
    currKIk = zeros(Complex{FT}, 3*myMolData.numatoms, 3*myMolData.numatoms);
    currKIblock = zeros(Complex{FT}, 3*myMolData.numatoms, 3*myMolData.numatoms);
    for nn2=1:myMolData.numblocks2
        for nn1=1:myMolData.numblocks1
            currKIblock[:, :] .= myMolData.KI[3*(nn1-1)*myMolData.numatoms+1:3*nn1*myMolData.numatoms,
                      3*(nn2-1)*myMolData.numatoms+1:3*nn2*myMolData.numatoms];
            lmul!(exp(-1im*dot(k,
                               myMolData.blocklist1[nn1] .*
                               myPeriodicData.latticevecs[:, 1] .+
                               myMolData.blocklist2[nn2] .*
                               myPeriodicData.latticevecs[:, 2])),
                  currKIblock);
            hermpart!(currKIblock);
            currKIk[:, :] .+= currKIblock[:, :];
        end
    end
    hermpart!(currKIk);
    return currKIk;
    
end

"""

    KIk!(myMolData::OneMol_WithPhonons{FT}, myPeriodicData::PeriodicData{FT},
         k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the `k`-dependence of the internuclear spring constant matrix for the molecule in
isolation from its real-space representation `myMolData.KI`, and store the result in
`myMolData.KIBloch`. The result will be the same as `myMolData.KI` if the molecule has no
periodic dimensions or couplings across unit cells.

See also: [`constructMolAlpha!`](@ref).

"""
function KIk!(myMolData::OneMol_WithPhonons{FT}, myPeriodicData::PeriodicData{FT},
              k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    fill!(myMolData.KIBloch, zero(eltype(myMolData.KIBloch)));
    if (myPeriodicData.numdims < 1)
        myMolData.KIBloch .= complex(myMolData.KI);
    else
        if (length(k) != 3)
            error("KIk(): k must have length 3");
        end
        currKIblock = zeros(Complex{FT}, 3*myMolData.numatoms, 3*myMolData.numatoms);
        for nn2=1:myMolData.numblocks2
            for nn1=1:myMolData.numblocks1
                currKIblock[:, :] .= myMolData.KI[3*(nn1-1)*myMolData.numatoms+1:3*nn1*myMolData.numatoms,
                                                  3*(nn2-1)*myMolData.numatoms+1:3*nn2*myMolData.numatoms];
                lmul!(exp(-1im*dot(k,
                                   myMolData.blocklist1[nn1] .*
                                   myPeriodicData.latticevecs[:, 1] .+
                                   myMolData.blocklist2[nn2] .*
                                   myPeriodicData.latticevecs[:, 2])),
                      currKIblock);
                hermpart!(currKIblock);
                myMolData.KIBloch[:, :] .+= currKIblock[:, :];
            end
        end
    end
    hermpart!(myMolData.KIBloch);
end
