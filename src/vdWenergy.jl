"""
    vdWenergy(myARGS)

Extract command-line arguments for computing the vdW interaction energy among a collection
of molecules, then compute the vdW energy integrand in imaginary frequency appropriately.

!!! note

Any frequencies will be converted to their absolute values times `1im`. Frequencies are not
allowed to be exactly zero.

!!! note

This computes the integrand ``\\ln(\\det(T_{\\infty} T^{-1}))`` at specified imaginary
frequencies, without doing any integration or Matsubara summation. This is typically a
negative dimensionless quantity, and no factors of ``2\\pi`` or ``\\hbar`` are included.

!!! tip

This function can be embarrassingly parallelized over frequency and wavevector. However, the
integrand at each frequency and wavevector is calculated for every transformation, as this
maximizes the efficiency of the program by virtue of reusing quantities independent of
geometric transformations.

!!! tip

The recommended way to do frequency integration or Matsubara summation is to supply a list
of frequencies (typically not too many), and then interpolate the integrand to do
integration or Matsubara summation after the fact.

"""
function vdWenergy(myARGS)

    mollistfilename = readRequiredFilenameARG(myARGS, "mollistfilename", "AARMBEMC_vdW.jl");
    outfilename = readRequiredFilenameARG(myARGS, "outfilename", "AARMBEMC_vdW.jl");
    
    translistfilename = readOptionalFilenameARG(myARGS, "translistfilename");
    periodicfilename = readOptionalFilenameARG(myARGS, "periodicfilename");
    
    Genvstr = readGenvstr(myARGS);
    
    nonretarded, outfilename = readOptionalBoolcondARG(myARGS, "nonretarded", outfilename);
    
    FT = readOptionalFloatType(myARGS);
    
    myAllMolData = readMolSystem(mollistfilename, false, FT);
    myTransData = readTransData(myAllMolData, translistfilename);
    myPeriodicData = readPeriodicData(FT, periodicfilename);

    myDeltaKISystemData = readDeltaKISystemData(mollistfilename, myAllMolData);
    
    freqlist = Array{FT, 1}(undef, 0);
    readfreqlist!(freqlist, myARGS, "AARMBEMC_vdW.jl");
    
    klist = readklist(myARGS, "AARMBEMC_vdW.jl", myPeriodicData.numdims, FT);

    freqklist = hcat(vec(repeat(transpose(freqlist), size(klist, 1), 1)),
                     repeat(klist, length(freqlist), 1));

    vdWenergy_general(outfilename, myAllMolData, myTransData, myPeriodicData,
                      myDeltaKISystemData, Genvstr, nonretarded, freqklist)
    
end

"""

    vdWenergy_general(outfilename, myAllMolData, myTransData, myPeriodicData,
                      myDeltaKISystemData, Genvstr, nonretarded, freqklist)

Compute the actual vdW energy integrands (for conduction, radiation, and both together) and
print each to `outfilename` (radiation integrand, then conduction integrand, then total
integrand).

"""
function vdWenergy_general(outfilename::AbstractString,
                           myAllMolData::MolSystem{<:OneMol{FT}},
                           myTransData::TransData{FT},
                           myPeriodicData::PeriodicData{FT},
                           myDeltaKISystemData::DeltaKISystemData{FT},
                           Genvstr::AbstractString, nonretarded::Bool,
                           freqklist::Array{FT, 2}) where FT<:AbstractFloat

    GFSCAGG! = GFVACGG_sca!;
    if (Genvstr == "PEC")
        GFSCAGG! = GFPECGG!;
    end

    totalYe0eeGenvmat = SharedArray{Complex{FT}, 2}(3*myAllMolData.cumnumatomslist[end],
                                                    3*myAllMolData.cumnumatomslist[end]);
    Ye0IIDeltaKImat = SharedArray{Complex{FT}, 2}(3*myAllMolData.cumnumatomslist[end],
                                                  3*myAllMolData.cumnumatomslist[end]);
    tempmat = SharedArray{Complex{FT}, 2}(3*myAllMolData.cumnumatomslist[end],
                                          3*myAllMolData.cumnumatomslist[end]);
    fill!(totalYe0eeGenvmat, zero(eltype(totalYe0eeGenvmat)));
    fill!(tempmat, zero(eltype(totalalphaeTinvmat)));
    fill!(Ye0IIDeltaKImat, zero(eltype(Ye0IIDeltaKImat)));

    radintegrandvaclist = zeros(FT, myAllMolData.nummols);
    condintegrandvaclist = zeros(FT, myAllMolData.nummols);
    radintegrandvac = zero(FT);
    condintegrandvac = zero(FT);
    radintegrandtot = zero(FT);
    condintegrandtot = zero(FT);
    bothintegrandtot = zero(FT);
  
    for ll=1:size(freqklist, 1)

        currfreq = 1im * abs(freqklist[ll, 1]);
        currfreqG = currfreq;
        if (nonretarded)
            currfreqG = zero(currfreqG);
        end
        currk = vec(freqklist[ll, 2:4]);
        fill!(totalYe0eeGenvmat, zero(eltype(totalYe0eeGenvmat)));
        fill!(tempmat, zero(eltype(totalalphaeTinvmat)));
        fill!(Ye0IIDeltaKImat, zero(eltype(Ye0IIDeltaKImat)));

        for nn=1:myAllMolData.nummols

            currumolind = myAllMolData.uniquemolinds[nn];
            if (currumolind > -1)
            
                constructMolZI0!(myAllMolData.uniqueMolDataArray[currumolind],
                                 myPeriodicData, currfreq, currk);

                constructMolYe0eeGvacinf!(myAllMolData.uniqueMolDataArray[currumolind],
                                          myPeriodicData, currfreqG, currk);
                                
                condintegrandvaclist[nn] =
                logabsdet(myAllMolData.uniqueMolDataArray[currumolind].KIbloch)[1];

                radintegrandvaclist[nn] =
                computeradintegrandvac(tempmat, myAllMolData.uniqueMolDataArray[currumolind]);

            end
            
        end

        radintegrandvac = integrandvac(myAllMolData, radintegrandvaclist);
        conddintegrandvac = integrandvac(myAllMolData, condintegrandvaclist);

        for tt=1:myTransData.numTrans

            fill!(tempmat, zero(eltype(totalalphaeTinvmat)));
            fill!(Ye0IIDeltaKImat, zero(eltype(Ye0IIDeltaKImat)));
          
            transformAtomPos!(myTransData, myAllMolData, tt);

            for nn2=1:myAllMolData.nummols
            
                startidx2 = 1;
                if (nn2 > 1)
                    startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
                end
                currmolind2 = myAllMolData.uniquemolinds[nn2];
                if (currmolind2 < 0)
                    currmolind2 = myAllMolData.dupemolinds[nn2];
                end

                if (myTransData.changedFromBefore[nn2, tt])
                    assembleYe0eeGenvDiagBlock!(totalYe0eeGenvmat, myAllMolData,
                                                myTransData, nn2, currmolind2, startidx2,
                                                tt, myPeriodicData, GFSCAGG!, currfreqG,
                                                currk);
                end

            end

            for nn2=1:myAllMolData.nummols
              
                startidx2 = 1;
                if (nn2 > 1)
                    startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
                end
        
                for nn1=(nn2+1):myAllMolData.nummols

                    startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;

                    if (myTransData.changedFromBefore[nn2, tt] ||
                        myTransData.changedFromBefore[nn1, tt])
                      
                        assembleYe0eeGenvOffDiagBlock!(totalYe0eeGenvmat, myAllMolData,
                                                       myTransData, nn1, nn2, startidx1,
                                                       startidx2, myPeriodicData, GFSCAGG!,
                                                       currfreqG, currk);
                    end
                end

            end

            for cc=1:length(myDeltaKISystemData.DeltaKImolblockdataarray)

                mol1 = myDeltaKISystemData.DeltaKImolblockdataarray[cc].mol1;
                mol2 = myDeltaKISystemData.DeltaKImolblockdataarray[cc].mol2;

                if (myTransData.changedFromBefore[mol2, tt] ||
                    myTransData.changedFromBefore[mol1, tt])

                    assembleDeltaKImolblock!(myDeltaKISystemData.DeltaKImolblockdataarray[cc],
                                             myAllMolData, myTransData, myPeriodicData);

                end

            end

            assembleDeltaKIdiagblocks!(myDeltaKISystemData, myAllMolData);

            assembleYe0IIDeltaKI!(Ye0IIDeltaKImat, myDeltaKISystemData, myAllMolData,
                                  myPeriodicData, currk);

            condintegrandtot = computecondintegrandwithDeltaKI(tempmat,
                                                               Ye0IIDeltaKImat,
                                                               myAllMolData,
                                                               myTransData,
                                                               tt) - condintegrandvac;

            radintegrandtot = computeradintegrandnoDeltaKI(tempmat,
                                                           totalYe0eeGenvmat,
                                                           myAllMolData,
                                                           myTransData,
                                                           tt) - radintegrandvac;

            bothintegrandtot = (computeradintegrandwithDeltaKI(tempmat,
                                                               totalYe0eeGenvmat,
                                                               Ye0IIDeltaKImat,
                                                               myAllMolData,
                                                               myTransData,
                                                               tt) - radintegrandvac) +
            condintegrandtot;

            outfile = open(outfilename, "a+");
            @printf(outfile, "%s %.16e ", myTransData.transLabels[tt], imag(currfreq));
            if (myPeriodicData.numdims > 0)
                @printf(outfile, "%.16e %.16e %.16e ", currk[1], currk[2], currk[3]);
            end
            @printf(outfile, "%.16e %.16e %.16e\n", radintegrandtot, condintegrandtot,
                    bothintegrandtot);
            if (myPeriodicData.numdims > 0)
                @printf(outfile, "%s %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
                        myTransData.transLabels[tt], imag(currfreq), -1*currk[1],
                        -1*currk[2], -1*currk[3], radintegrandtot, condintegrandtot,
                        bothintegrandtot);
            end
            close(outfile);
            
            untransformAtomPos!(myTransData, myAllMolData, tt);

        end
    
    end
  
end

"""

    constructMolZI0!(myMolData::OneMol_WithPhonons{FT},
                     myPeriodicData::PeriodicData{FT},
                     freq::Union{FT, Complex{FT}},
                     k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the product of the uncoupled nuclear response with the internuclear spring
constant matrix of a single molecule. This also constructs the atomic polarizabilities used
to evaluate the Gaussian widths used for the Green's function matrices.

"""
function constructMolZI0!(myMolData::OneMol_WithPhonons{FT},
                          myPeriodicData::PeriodicData{FT},
                          freq::Union{FT, Complex{FT}},
                          k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    constructMolAlpha!(myMolData, myPeriodicData, freq, k);
    mulMatDiagleft!(myMolData.KIbloch, vecNto3N(myMolData.Ye0II));
    addtodiag!(myMolData.KIbloch, ones(eltype(myMolData.KIbloch), 3*myMolData.numatoms));

end

"""

    constructMolYe0eeGvacinf!(myMolData::OneMol_WithPhonons{FT},
                              myPeriodicData::PeriodicData{FT},
                              currfreqG::Union{FT, Complex{FT}},
                              currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the Green's function interaction matrix of a molecule in vacuum
`myMolData.GvacGGinf`, then multiplies on the left by its electronic response in-place, at
frequency `currfreqG` and wavevector `currk` passing `myMolData` and `myPeriodicData` to the
subsidiary functions. This matrix is then subtracted from the identity.

"""
function constructMolYe0eeGvacinf!(myMolData::OneMol_WithPhonons{FT},
                                   myPeriodicData::PeriodicData{FT},
                                   currfreqG::Union{FT, Complex{FT}},
                                   currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat
  
    assemble1MolGvacinf!(myMolData, myPeriodicData, currfreqG, currk);
    mulMatDiagleft!(myMolData.GvacGGinf, vecNto3N(myMolData.Ye0ee));
    subtractfromdiag!(myMolData.GvacGGinf, ones(eltype(myMolData.GvacGGinf),
                                                3*myMolData.numatoms));

end

"""

    computeradintegrandvac(tempmat::SharedArray{Complex{FT}, 2},
                           myMolData::OneMol_WithPhonons{FT}) where FT<:AbstractFloat

Set up the matrix `tempmat` to compute the vacuum radiative contribution to the vdW energy
integrand.

"""
function computeradintegrandvac(tempmat::SharedArray{Complex{FT}, 2},
                                myMolData::OneMol_WithPhonons{FT}) where FT<:AbstractFloat
  
    fill!(tempmat, zero(eltype(tempmat)));
                
    tempmat[1:3*myMolData.numatoms, 1:3*myMolData.numatoms] .= inv(myMolData.KIbloch);

    tempview = view(tempmat, 1:3*myMolData.numatoms, 1:3*myMolData.numatoms);
    
    lmul!(tempview, -1);
    mulMatDiagright!(tempview, vecNto3N(myMolData.Ye0II .* myMolData.Ke));
    mulMatDiagleft!(tempview, vecNto3N(myMolData.Ye0ee .* myMolData.Ke));
    tempview .+= myMolData.GvacGGinf;

end
    
"""

    integrandvac(myAllMolData::MolSystem{<:OneMol{FT}},
                 integrandvaclist::Array{FT, 1}) where FT<:AbstractFloat

Populate the list of vacuum contributions to the vdW integrand `integrandvaclist` with
contributions from duplicate molecules where they exist, then return the sum yielding the
total contribution to the integrand in vacuum.

"""
function integrandvac(myAllMolData::MolSystem{<:OneMol{FT}},
                      integrandvaclist::Array{FT, 1}) where FT<:AbstractFloat

    myintegrandvac = zero(eltype(integrandvaclist));
    for nn=1:myAllMolData.nummols
        if (myAllMolData.dupemolinds[nn] > 0)
            integrandvaclist[nn] = integrandvaclist[myAllMolData.dupemolinds[nn]];
        end
        myintegrandvac += integrandvaclist[nn];
    end

    return myintegrandvac;

end

"""

    assembleYe0eeGenvDiagBlock!(totalYe0eeGenvmat::SharedArray{Complex{FT}, 2},
                                myAllMolData::MolSystem{<:OneMol{FT}},
                                myTransData::TransData{FT}, nn2::Integer,
                                currmolind2::Integer, startidx2::Integer,
                                tt::Integer, myPeriodicData::PeriodicData{FT},
                                GFSCAGG!::Function,
                                currfreqG::Union{FT, Complex{FT}},
                                currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the molecular diagonal block of the scattering Green's function interaction matrix
of molecule `nn2` (whose unique molecular index `currmolind2` is precomputed and passed as
an argument) at frequency `currfreqG` and wavevector `currk` using `GFSCAGG!`, and multiply
on the left by the electronic response, then add the vacuum contribution which has already
been computed as the electronic polarizabilities multiplied by the vacuum Green's function
matrix of molecule `nn2`, performing any required matrix rotations corresponding to
transformation `tt` for molecule `nn2`, and stamp that into
`totalYe0eeGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
Finally, undo any rotation corresponding to transformation `tt` in the vacuum Green's
function matrix of molecule `nn2` for reuse with future transformations.

"""
function assembleYe0eeGenvDiagBlock!(totalYe0eeGenvmat::SharedArray{Complex{FT}, 2},
                                     myAllMolData::MolSystem{<:OneMol{FT}},
                                     myTransData::TransData{FT}, nn2::Integer,
                                     currmolind2::Integer, startidx2::Integer,
                                     tt::Integer, myPeriodicData::PeriodicData{FT},
                                     GFSCAGG!::Function,
                                     currfreqG::Union{FT, Complex{FT}},
                                     currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalYe0eeGenvmat, (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));

    fill!(tempview, zero(eltype(tempview)));

    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn2, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);

    lmul!(-1, tempview);
    mulMatDiagleft!(tempview, vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].Ye0ee));

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                 myTransData.baselineRot[3*nn2-2:3*nn2, :]);
    
    tempview .+= myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf[:, :];

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 transpose(myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                           myTransData.baselineRot[3*nn2-2:3*nn2, :]));

end

"""

    assembleYe0eeGenvOffDiagBlock!(totalYe0eeGenvmat::SharedArray{Complex{FT}, 2},
                                   myAllMolData::MolSystem{<:OneMol{FT}},
                                   myTransData::TransData{FT},
                                   nn1::Integer, nn2::Integer,
                                   startidx1::Integer,
                                   startidx2::Integer,
                                   myPeriodicData::PeriodicData{FT},
                                   GFSCAGG!::Function,
                                   currfreqG::Union{FT, Complex{FT}},
                                   currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the molecular off-diagonal block of the scattering Green's function interaction
matrix of molecules `nn2` as the source and `nn1` as the field (whose respective unique
molecular indices `currmolind1` and `currmolind2` are precomputed and passed as arguments)
at frequency `currfreqG` and wavevector `currk` using `GFSCAGG!`, multiply on
the left by the electronic response, add the vacuum contribution, and stamp that into
`totalYe0eeGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
Do the same for the adjoint, swapping `nn2` and `nn1` without performing further Green's
function computations.

"""
function assembleYe0eeGenvOffDiagBlock!(totalYe0eeGenvmat::SharedArray{Complex{FT}, 2},
                                        myAllMolData::MolSystem{<:OneMol{FT}},
                                        myTransData::TransData{FT},
                                        nn1::Integer, nn2::Integer,
                                        startidx1::Integer,
                                        startidx2::Integer,
                                        myPeriodicData::PeriodicData{FT},
                                        GFSCAGG!::Function,
                                        currfreqG::Union{FT, Complex{FT}},
                                        currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalYe0eeGenvmat, (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));
    fill!(tempview, zero(eltype(tempview)));

    assembleGvacMolBlock!(tempview, myAllMolData, myTransData, nn1, nn2, 1, 1,
                          myPeriodicData, currfreqG, currk);
    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn1, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);
    lmul!(-1, tempview);

    totalYe0eeGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                      (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1])] .=
    adjoint(totalYe0eeGenvmat[(startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]),
                              (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])])

    currmolind2 = myAllMolData.uniquemolinds[nn2];
    if (currmolind2 < 0)
        currmolind2 = myAllMolData.dupemolinds[nn2];
    end

    currmolind1 = myAllMolData.uniquemolinds[nn1];
    if (currmolind1 < 0)
        currmolind1 = myAllMolData.dupemolinds[nn1];
    end

    mulMatDiagleft!(tempview, vecNto3N(myAllMolData.uniqueMolDataArray[currmolind1].Ye0ee));
    
    tempview = view(totalYe0eeGenvmat, (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                    (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]));

    mulMatDiagleft!(tempview, vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].Ye0ee));
    
end

"""

    assembleYe0IIDeltaKI!(Ye0IIDeltaKImat::AbstractArray{Complex{FT}, 2},
                          myDeltaKISystemData::DeltaKISystemData{FT},
                          myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                          blochk::AbstractArray{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the product of the uncoupled nuclear response across all molecules with the
internuclear spring constant matrix coupling different molecules.

"""
function assembleYe0IIDeltaKI!(Ye0IIDeltaKImat::AbstractArray{Complex{FT}, 2},
                               myDeltaKISystemData::DeltaKISystemData{FT},
                               myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                               blochk::AbstractArray{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    assembleDeltaKI!(Ye0IIDeltaKImat, myDeltaKISystemData, myAllMolData, currk);
    for nn=1:myAllMolData.nummols

        startidx = 1;
        if (nn > 1)
            startidx = 3*myAllMolData.cumnumatomslist[nn - 1] + 1;
        end
        tempview = view(Ye0IIDeltaKImat,
                        startidx:3*myAllMolData.cumnumatomslist[nn],
                        1:3*myAllMolData.cumnumatomslist[end]);

        currumolind = myAllMolData.uniquemolinds[nn];
        if (currumolind < 0)
            currumolind = myAllMolData.dupemolinds[nn]
        end

        mulMatDiagleft!(tempview,
                        vecNto3N(myAllMolData.uniqueMolDataArray[currumolind].Ye0II));

    end
    
end

"""

    computeradintegrandnoDeltaKI(tempmat::SharedArray{Complex{FT}, 2},
                                 totalYe0eeGenvmat::SharedArray{Complex{FT}, 2},
                                 myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                 myTransData::TransData{FT},
                                 tt::Integer) where FT<:AbstractFloat

Compute and return the total contribution (before subtracting the vacuum contribution) to
the vdW energy integrand from radiation in the absence of conductive couplings between
different molecules, at transformation index `tt`.

"""
function computeradintegrandnoDeltaKI(tempmat::SharedArray{Complex{FT}, 2},
                                      totalYe0eeGenvmat::SharedArray{Complex{FT}, 2},
                                      myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                      myTransData::TransData{FT},
                                      tt::Integer) where FT<:AbstractFloat

    fill!(tempmat, zero(eltype(tempmat)));

    for nn=1:myAllMolData.nummols

        currumolind = myAllMolData.uniquemolinds[nn];
        if (currumolind < 0)
            currumolind = myAllMolData.dupemolinds[nn];
        end

        rotate3Nmat!(myAllMolData.uniqueMolDataArray[currumolind].KIBloch,
                     myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] *
                     myTransData.baselineRot[3*nn-2:3*nn, :]);

        startidx = 1;
        if (nn > 1)
            startidx = 3*myAllMolData.cumnumatomslist[nn-1] + 1;
        end
        endidx = 3*myAllMolData.cumnumatomslist[nn];

        tempview = view(tempmat, startidx:endidx, startidx:endidx);
        fill!(tempview, zero(eltype(tempview)));
        tempview .+= inv(myAllMolData.uniqueMolDataArray[currumolind].KIBloch);
        lmul!(tempview, -1);
        mulMatDiagright!(tempview,
                         vecNto3N(myAllMolData.uniqueMolDataArray[currumolind].Ye0II .*
                                  myAllMolData.uniqueMolDataArray[currumolind].Ke));
        mulMatDiagleft!(tempview,
                        vecNto3N(myAllMolData.uniqueMolDataArray[currumolind].Ye0ee .*
                                 myAllMolData.uniqueMolDataArray[currumolind].Ke));
        
        rotate3Nmat!(myAllMolData.uniqueMolDataArray[currumolind].KIBloch,
                     transpose(myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] *
                               myTransData.baselineRot[3*nn-2:3*nn, :]));

    end

    tempmat .+= totalYe0eeGenvmat;
    return logabsdet(tempmat)[1];

end

"""

    computeradintegrandwithDeltaKI(tempmat::SharedArray{Complex{FT}, 2},
                                   totalYe0eeGenvmat::SharedArray{Complex{FT}, 2},
                                   Ye0IIDeltaKImat::SharedArray{Complex{FT}, 2},
                                   myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                   myTransData::TransData{FT},
                                   tt::Integer) where FT<:AbstractFloat

Compute and return the total contribution (before subtracting the vacuum contribution) to
the vdW energy integrand from radiation in the presence of conductive couplings between
different molecules, at transformation index `tt`.

"""
function computeradintegrandwithDeltaKI(tempmat::SharedArray{Complex{FT}, 2},
                                        totalYe0eeGenvmat::SharedArray{Complex{FT}, 2},
                                        Ye0IIDeltaKImat::SharedArray{Complex{FT}, 2},
                                        myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                        myTransData::TransData{FT},
                                        tt::Integer) where FT<:AbstractFloat

    fill!(tempmat, zero(eltype(tempmat)));

    for nn=1:myAllMolData.nummols

        currumolind = myAllMolData.uniquemolinds[nn];
        if (currumolind < 0)
            currumolind = myAllMolData.dupemolinds[nn];
        end

        rotate3Nmat!(myAllMolData.uniqueMolDataArray[currumolind].KIBloch,
                     myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] *
                     myTransData.baselineRot[3*nn-2:3*nn, :]);

        startidx = 1;
        if (nn > 1)
            startidx = 3*myAllMolData.cumnumatomslist[nn-1] + 1;
        end
        endidx = 3*myAllMolData.cumnumatomslist[nn];

        tempview = view(tempmat, startidx:endidx, startidx:endidx);
        fill!(tempview, zero(eltype(tempview)));
        tempview .+= myAllMolData.uniqueMolDataArray[currumolind].KIBloch;
        
        rotate3Nmat!(myAllMolData.uniqueMolDataArray[currumolind].KIBloch,
                     transpose(myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] *
                               myTransData.baselineRot[3*nn-2:3*nn, :]));

    end

    tempmat .+= Ye0IIDeltaKImat;
    tempmat .= inv(tempmat);
    lmul!(tempmat, -1);

    for nn2=1:myAllMolData.nummols

        currumolind2 = myAllMolData.uniquemolinds[nn2];
        if (currumolind2 < 0)
            currumolind2 = myAllMolData.dupemolinds[nn2];
        end

        startidx2 = 1;
        if (nn2 > 1)
            startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
        end
        endidx2 = 3*myAllMolData.cumnumatomslist[nn2];

        for nn1=1:myAllMolData.nummols

            currumolind1 = myAllMolData.uniquemolinds[nn1];
            if (currumolind1 < 0)
                currumolind1 = myAllMolData.dupemolinds[nn1];
            end
            
            startidx1 = 1;
            if (nn1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
            end
            endidx1 = 3*myAllMolData.cumnumatomslist[nn1];

            tempview = view(tempmat, startidx1:endidx1, startidx2:endidx2);

            mulMatDiagright!(tempview,
                             vecNto3N(myAllMolData.uniqueMolDataArray[currumolind2].Ye0II .*
                                      myAllMolData.uniqueMolDataArray[currumolind2].Ke));
            mulMatDiagleft!(tempview,
                            vecNto3N(myAllMolData.uniqueMolDataArray[currumolind1].Ye0ee .*
                                     myAllMolData.uniqueMolDataArray[currumolind1].Ke));
            
        end
    end

    tempmat .+= totalYe0eeGenvmat;
    return logabsdet(tempmat)[1];

end

"""

    computecondintegrandwithDeltaKI(tempmat::SharedArray{Complex{FT}, 2},
                                    Ye0IIDeltaKImat::SharedArray{Complex{FT}, 2},
                                    myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                    myTransData::TransData{FT},
                                    tt::Integer) where FT<:AbstractFloat

Compute and return the total contribution (with the vacuum contribution vanishing) to the
vdW energy integrand from conduction (in the presence of conductive couplings between
different molecules), at transformation index `tt`.

"""
function computecondintegrandwithDeltaKI(tempmat::SharedArray{Complex{FT}, 2},
                                         Ye0IIDeltaKImat::SharedArray{Complex{FT}, 2},
                                         myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                         myTransData::TransData{FT},
                                         tt::Integer) where FT<:AbstractFloat

    fill!(tempmat, zero(eltype(tempmat)));

    for nn=1:myAllMolData.nummols

        currumolind = myAllMolData.uniquemolinds[nn];
        if (currumolind < 0)
            currumolind = myAllMolData.dupemolinds[nn];
        end

        rotate3Nmat!(myAllMolData.uniqueMolDataArray[currumolind].KIBloch,
                     myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] *
                     myTransData.baselineRot[3*nn-2:3*nn, :]);

        startidx = 1;
        if (nn > 1)
            startidx = 3*myAllMolData.cumnumatomslist[nn-1] + 1;
        end
        endidx = 3*myAllMolData.cumnumatomslist[nn];

        tempview = view(tempmat, startidx:endidx, startidx:endidx);
        fill!(tempview, zero(eltype(tempview)));
        tempview .+= myAllMolData.uniqueMolDataArray[currumolind].KIBloch;
        
        rotate3Nmat!(myAllMolData.uniqueMolDataArray[currumolind].KIBloch,
                     transpose(myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] *
                               myTransData.baselineRot[3*nn-2:3*nn, :]));

    end

    tempmat .+= Ye0IIDeltaKImat;
    return logabsdet(tempmat)[1];

end


## TO DELETE: constructMolalphaealphainv!(...) because the relevant
## quantities for combined conduction and radiation are different
"""

    constructMolalphaealphainv!(myMolData::OneMol_WithPhonons{FT},
                                myPeriodicData::PeriodicData{FT},
                                freq::Union{FT, Complex{FT}},
                                k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the product of the electronic susceptibility with the full inverse molecular
susceptibility `myMolData.alphainv` in-place at frequency `freq` and wavevector `k`, the
latter of which must be a 3-element real vector. With phonons, `myMolData.alphainv` will
be a matrix, the vector of atomic polarizabilities `myMolData.alpha0` will depend on the
elements of `myMolData.alphainv`, and `myMolData.alphae` will be independent of these two.
Without phonons, all three will be simply related to each other. This is more efficient
than calling `constructMolAlpha!` and then performing a `broadcast` multiplication involving
`myMolData.alphae` multiplying the rows of `myMolData.alphainv`.

"""
function constructMolalphaealphainv!(myMolData::OneMol_NoPhonons{FT},
                                     myPeriodicData::PeriodicData{FT},
                                     freq::Union{FT, Complex{FT}},
                                     k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    setindex!(myMolData.alphae,
              (myMolData.Qe.^2) ./ (myMolData.Ke .- (1im.*freq.*myMolData.Be .+
                                freq^2 .* myMolData.Me)),
              1:myMolData.numatoms);
    setindex!(myMolData.alpha0, abs.(myMolData.alphae), 1:myMolData.numatoms);
    fill!(myMolData.alphainv, one(eltype(myMolData.alphainv)));

end
function constructMolalphaealphainv!(myMolData::OneMol_WithPhonons{FT},
                                     myPeriodicData::PeriodicData{FT},
                                     freq::Union{FT, Complex{FT}},
                                     k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat
  
    Ke3N = vecNto3N(myMolData.Ke);
    Qe3N = vecNto3N(myMolData.Qe);
  
    setindex!(myMolData.alphainv, KIk(myMolData, myPeriodicData, k),
              1:3*myMolData.numatoms, 1:3*myMolData.numatoms);
    
    addtodiag!(myMolData.alphainv,
               Ke3N .- vecNto3N(1im.*freq.*myMolData.BI .+ freq^2 .* myMolData.MI));

    setindex!(myMolData.alphainv, inv(myMolData.alphainv),
              1:3*myMolData.numatoms, 1:3*myMolData.numatoms);

    mulMatDiagleft!(myMolData.alphainv, Ke3N);
    mulMatDiagright!(myMolData.alphainv, Ke3N);

    alpha = deepcopy(myMolData.alphainv);
    subtractfromdiag!(alpha,
                      Ke3N .- vecNto3N(1im.*freq.*myMolData.Be .+ freq^2 .* myMolData.Me));

    mulMatDiagright!(alpha, 1 ./ Qe3N);
    mulMatDiagleft!(alpha, 1 ./ Qe3N);

    alpha = inv(alpha);
    alphacarttr = Array{eltype(alpha), 2}(undef, myMolData.numatoms, myMolData.numatoms);
    for qq=1:myMolData.numatoms
        for pp=1:myMolData.numatoms
            alphacarttr[pp,qq] = tr(alpha[3*pp-2:3*pp, 3*qq-2:3*qq]);
        end
    end
    setindex!(myMolData.alpha0, vec(abs.(sum(alphacarttr, dims=2))), 1:myMolData.numatoms);
    lmul!(1/3, myMolData.alpha0);

    mulMatDiagright!(myMolData.alphainv, 1 ./ Qe3N);
    
    mulMatDiagleft!(myMolData.alphainv,
                    Qe3N  ./
                    (Ke3N .- vecNto3N(1im.*freq.*myMolData.Be .+
                                      freq^2 .* myMolData.Me)));
    
    subtractfromdiag!(myMolData.alphainv, ones(eltype(myMolData.alphainv),
                                               3*myMolData.numatoms));

    setindex!(myMolData.alphae,
              myMolData.Qe.^2 ./ (myMolData.Ke .-
                                  (1im.*freq.*myMolData.Be .+
                                   freq^2 .* myMolData.Me)),
              1:myMolData.numatoms);
  
end

## TO DELETE: constructMolalphaeGvacinf!(...) because the nomenclature
## should reflect the variables used
"""

    constructMolalphaeGvacinf!(myMolData::OneMol{FT},
                               myPeriodicData::PeriodicData{FT},
                               currfreqG::Union{FT, Complex{FT}},
                               currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the Green's function interaction matrix of a molecule in vacuum
`myMolData.GvacGGinf`, then multiplies on the left by its electronic polarizabilities
in-place, at frequency `currfreqG` and wavevector `currk` passing `myMolData` and
`myPeriodicData` to the subsidiary functions.

"""
function constructMolalphaeGvacinf!(myMolData::OneMol{FT},
                                    myPeriodicData::PeriodicData{FT},
                                    currfreqG::Union{FT, Complex{FT}},
                                    currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat
  
    assemble1MolGvacinf!(myMolData, myPeriodicData, currfreqG, currk);
    mulMatDiagleft!(myMolData.GvacGGinf, vecNto3N(myMolData.alphae));

end

## TO DELETE: mulalphaeGBlock!(...) because the relevant quantities
## have changed
"""

    mulalphaeGBlock!(A::AbstractArray{Complex{FT}, 2},
                     myAllMolData::MolSystem{<:OneMol{FT}},
                     nn1::Integer, startidx1::Integer,
                     startidx2::Integer,
                     endidx2::Integer) where FT<:AbstractFloat

Assuming `A` represents a Green's function interaction matrix, multiply its block
`A[startidx1-1+3*myAllMolData.numatomslist[nn1], startidx2:endidx2]` on the left
by the electronic polarizabilities of molecule `nn1` from `myAllMolData`.

"""
function mulalphaeGBlock!(A::AbstractArray{Complex{FT}, 2},
                          myAllMolData::MolSystem{<:OneMol{FT}},
                          nn1::Integer, startidx1::Integer,
                          startidx2::Integer,
                          endidx2::Integer) where FT<:AbstractFloat

    currmolind1 = myAllMolData.dupemolinds[nn1];
    if (currmolind1 < 0)
        currmolind1 = myAllMolData.uniquemolinds[nn1];
    end

    if (size(A, 1) < 3*myAllMolData.numatomslist[nn1])
        error("A (currently with ", size(A, 1), " rows) should have at least ",
              3*myAllMolData.numatomslist[nn1], " rows in mulalphaeGBlock!()");
    end

    if (startidx1 < 1 || startidx1 > size(A, 1) - 3*myAllMolData.numatomslist[nn1] + 1)
        error("mulalphaeGBlock!(): startidx1 (currently ", startidx1,
              ") must be between 1 and ", size(A, 1) - 3*myAllMolData.numatomslist[nn1] + 1,
              " as A has ", size(A, 1), " rows and the block will have ",
              3*myAllMolData.numatomslist[nn1], " rows");
    end
    if (startidx2 < 1 || startidx2 > size(A, 2))
        error("mulalphaeGBlock!(): startidx2 (currently ", startidx2,
              ") must be between 1 and ", size(A, 2),
              " as A has ", size(A, 2), " columns");
    end
    if (endidx2 < 1 || endidx2 > size(A, 2))
        error("mulalphaeGBlock!(): endidx2 (currently ", startidx2,
              ") must be between 1 and ", size(A, 2),
              " as A has ", size(A, 2), " columns");
    end
    if (startidx2 > endidx2)
        error("mulalphaeGBlock!(): startidx2 (currently ", startidx2,
              ") must not be larger than endidx2 (currently ", endidx2, ")");
    end

    tempview = view(A, (startidx1-1).+1:3*myAllMolData.numatomslist[nn1], startidx2:endidx2);
    mulMatDiagleft!(tempview, vecNto3N(myAllMolData.uniqueMolDataArray[currmolind1].alphae));

end

## TO DELETE: assemblealphaeTinvDiagBlock!(...) because the relevant
## quantities have changed
"""

    assemblealphaeTinvDiagBlock!(totalalphaeTinvmat::SharedArray{Complex{FT}, 2},
                                 myAllMolData::MolSystem{<:OneMol{FT}},
                                 myTransData::TransData{FT}, nn2::Integer,
                                 currmolind2::Integer, startidx2::Integer,
                                 tt::Integer, myPeriodicData::PeriodicData{FT},
                                 GFSCAGG!::Function,
                                 currfreqG::Union{FT, Complex{FT}},
                                 currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the molecular diagonal block of the scattering Green's function interaction matrix
of molecule `nn2` (whose unique molecular index `currmolind2` is precomputed and passed as
an argument) at frequency `currfreqG` and wavevector `currk` using `GFSCAGG!`, multiply on
the left by the electronic polarizabilities, then add the vacuum contribution which has
already been computed as the electronic polarizabilities multiplied by the inverse
T-operator of molecule `nn2`, performing any required matrix rotations corresponding to
transformation `tt` for molecule `nn2`, and stamp that into
`totalalphaeTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
Finally, undo any rotation corresponding to transformation `tt` in the inverse T-operator of
molecule `nn2` for reuse with future transformations.

"""
function assemblealphaeTinvDiagBlock!(totalalphaeTinvmat::SharedArray{Complex{FT}, 2},
                                      myAllMolData::MolSystem{<:OneMol{FT}},
                                      myTransData::TransData{FT}, nn2::Integer,
                                      currmolind2::Integer, startidx2::Integer,
                                      tt::Integer, myPeriodicData::PeriodicData{FT},
                                      GFSCAGG!::Function,
                                      currfreqG::Union{FT, Complex{FT}},
                                      currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalalphaeTinvmat, (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));

    fill!(tempview, zero(eltype(tempview)));

    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn2, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);

    lmul!(-1, tempview);
    mulalphaeGBlock!(tempview, myAllMolData, nn2, 1, 1, 3*myAllMolData.numatomslist[nn2]);

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                 myTransData.baselineRot[3*nn2-2:3*nn2, :]);
    
    tempview .+= myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf[:, :];

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 transpose(myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                           myTransData.baselineRot[3*nn2-2:3*nn2, :]));

end

## TO DELETE: assemblealphaeTinvOffDiagBlock!(...) because the
## relevant quantities have changed
"""

    assemblealphaeTinvOffDiagBlock!(totalalphaeTinvmat::SharedArray{Complex{FT}, 2},
                                    myAllMolData::MolSystem{<:OneMol{FT}},
                                    myTransData::TransData{FT},
                                    nn1::Integer, nn2::Integer,
                                    startidx1::Integer,
                                    startidx2::Integer,
                                    myPeriodicData::PeriodicData{FT},
                                    GFSCAGG!::Function,
                                    currfreqG::Union{FT, Complex{FT}},
                                    currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the molecular off-diagonal block of the scattering Green's function interaction
matrix of molecules `nn2` as the source and `nn1` as the field (whose respective unique
molecular indices `currmolind1` and `currmolind2` are precomputed and passed as arguments)
at frequency `currfreqG` and wavevector `currk` using `GFSCAGG!`, multiply on
the left by the electronic polarizabilities, add the vacuum contribution, and stamp the
negative of that into
`totalalphaeTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
Do the same for the adjoint, swapping `nn2` and `nn1` without performing further Green's
function computations.

"""
function assemblealphaeTinvOffDiagBlock!(totalalphaeTinvmat::SharedArray{Complex{FT}, 2},
                                         myAllMolData::MolSystem{<:OneMol{FT}},
                                         myTransData::TransData{FT},
                                         nn1::Integer, nn2::Integer,
                                         startidx1::Integer,
                                         startidx2::Integer,
                                         myPeriodicData::PeriodicData{FT},
                                         GFSCAGG!::Function,
                                         currfreqG::Union{FT, Complex{FT}},
                                         currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalalphaeTinvmat, (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));
    fill!(tempview, zero(eltype(tempview)));

    assembleGvacMolBlock!(tempview, myAllMolData, myTransData, nn1, nn2, 1, 1,
                          myPeriodicData, currfreqG, currk);
    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn1, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);
    lmul!(-1, tempview);

    totalalphaeTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                       (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1])] .=
    adjoint(totalalphaeTinvmat[(startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]),
                               (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])])

    
    ## setindex!(totalalphaeTinvmat, adjoint(tempview),
    ##           (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
    ##           (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]));

    ## mulalphaeGBlock!(totalalphaeTinvmat, myAllMolData, nn1, startidx1, startidx2,
    ##                  startidx2 - 1 + 3*myAllMolData.numatomslist[nn2]);
    ## mulalphaeGBlock!(totalalphaeTinvmat, myAllMolData, nn2, startidx2, startidx1,
    ##                  startidx1 - 1 + 3*myAllMolData.numatomslist[nn1]);

    mulalphaeGBlock!(tempview, myAllMolData, nn1, 1, 1, 3*myAllMolData.numatomslist[nn2]);
    tempview = view(totalalphaeTinvmat, (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                    (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]));
    mulalphaeGBlock!(tempview, myAllMolData, nn2, 1, 1, 3*myAllMolData.numatomslist[nn1]);

end

