"""

    heattransfer(myARGS)

Extract command-line arguments for computing the radiative heat transfer among a collection
of molecules, then compute the radiative heat transfer appropriately.

!!! note

Any frequencies will be converted to their absolute values. Frequencies are not allowed to
be exactly zero.

!!! note

This computes the integrand ``\\Phi^{(m)}_{n}`` at specified real frequencies, without doing
any integration. This includes neither numerical prefactors outside of the trace, nor the
Planck function.

!!! tip

This function can be embarrassingly parallelized over frequency and wavevector. However, the
integrand at each frequency and wavevector is calculated for every transformation, as this
maximizes the efficiency of the program by virtue of reusing quantities independent of
geometric transformations.

"""
function heattransfer(myARGS)

    mollistfilename = readRequiredFilenameARG(myARGS, "mollistfilename",
                                              "AARMBEMC_heat.jl");
    outfilename = readRequiredFilenameARG(myARGS, "outfilename", "AARMBEMC_heat.jl");
  
    translistfilename = readOptionalFilenameARG(myARGS, "translistfilename");
    periodicfilename = readOptionalFilenameARG(myARGS, "periodicfilename");

    Genvstr = readGenvstr(myARGS);
    
    nonretarded, outfilename = readOptionalBoolcondARG(myARGS, "nonretarded", outfilename);

    FT = readOptionalFloatType(myARGS);
    println("FT is ", FT);
    
    myAllMolData = readMolSystem(mollistfilename, false, FT);
    myTransData = readTransData(myAllMolData, translistfilename);
    myPeriodicData = readPeriodicData(FT, periodicfilename);

    myDeltaKISystemData = readDeltaKISystemData(mollistfilename, myAllMolData);
    
    freqlist = Array{FT, 1}(undef, 0);
    readfreqlist!(freqlist, myARGS, "AARMBEMC_heat.jl");

    klist = readklist(myARGS, "AARMBEMC_heat.jl", myPeriodicData.numdims, FT);
    
    freqklist = hcat(vec(repeat(transpose(freqlist), size(klist, 1), 1)),
                     repeat(klist, length(freqlist), 1));

    heattransfer_general(outfilename, myAllMolData, myTransData, myPeriodicData,
                         myDeltaKISystemData, Genvstr, nonretarded, freqklist)
    
end

"""

    heattransfer_general(outfilename, myAllMolData, myTransData, myPeriodicData,
                         myDeltaKISystemData, Genvstr, nonretarded, freqklist)

Compute the actual heat transfer integrands (for conduction, radiation, and both together)
and print each to appropriate files based on the name specified by `outfilename`.

"""
function heattransfer_general(outfilename::AbstractString,
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

    totalGenvmat = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                            3*myAllMolData.cumnumatomslist[end]);
    totalGenvmatwithdiag = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                                    3*myAllMolData.cumnumatomslist[end]);
    DeltaKImat = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                          3*myAllMolData.cumnumatomslist[end]);
    totalTeemat = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                           3*myAllMolData.cumnumatomslist[end]);
    totalZIIinvmat = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                              3*myAllMolData.cumnumatomslist[end]);
    tempmat = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                              3*myAllMolData.cumnumatomslist[end]);
    Phimn = zeros(FT, myAllMolData.nummols, myAllMolData.nummols);
    
    for ll=1:size(freqklist, 1)

        currfreq = abs(freqklist[ll, 1]);

        currfreqG = currfreq;
        if (nonretarded)
            currfreqG = zero(currfreqG);
        end
        currk = vec(freqklist[ll, 2:4]);
        fill!(totalGenvmat, zero(eltype(totalGenvmat)));
        fill!(totalGenvmatwithdiag, zero(eltype(totalGenvmatwithdiag)));
        fill!(totalTeemat, zero(eltype(totalTeemat)));
        fill!(DeltaKImat, zero(eltype(DeltaKImat)));
        fill!(totalZIIinvmat, zero(eltype(totalZIIinvmat)));
        fill!(tempmat, zero(eltype(tempmat)));

        for nn=1:myAllMolData.nummols

            currumolind = myAllMolData.uniquemolinds[nn];

            if (currumolind > -1)
            
                constructMolAlpha!(myAllMolData.uniqueMolDataArray[currumolind],
                                   myPeriodicData, currfreq, currk);
                
                assemble1MolGvacinf!(myAllMolData.uniqueMolDataArray[currumolind],
                                     myPeriodicData, currfreqG, currk);

                assemble1MolGvacinfdiag!(myAllMolData.uniqueMolDataArray[currumolind],
                                         currfreqG);                
                
            end
            
        end

        for tt=1:myTransData.numTrans

            fill!(totalTeemat, zero(eltype(totalTeemat)));
            fill!(DeltaKImat, zero(eltype(DeltaKImat)));
            fill!(totalZIIinvmat, zero(eltype(totalZIIinvmat)));
            fill!(tempmat, zero(eltype(tempmat)));
          
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
                    assembleGenvDiagBlock!(totalGenvmat, totalGenvmatwithdiag, myAllMolData,
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

                        assembleGenvOffDiagBlock!(totalGenvmat, totalGenvmatwithdiag,
                                                  myAllMolData, myTransData, nn1, nn2,
                                                  startidx1, startidx2, myPeriodicData,
                                                  GFSCAGG!, currfreqG, currk);
                    end
                end

            end

            if (myPeriodicData.numdims > 0)
            
                for nn2=1:myAllMolData.nummols
                  
                    startidx2 = 1;
                    if (nn2 > 1)
                        startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
                    end

                    for nn1=(nn2+1):myAllMolData.nummols
                      
                        startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
                        
                        if (myTransData.changedFromBefore[nn2, tt] ||
                            myTransData.changedFromBefore[nn1, tt])
                          
                            assembleGenvOffDiagBlock!(totalGenvmat, totalGenvmatwithdiag,
                                                      myAllMolData, myTransData, nn2, nn1,
                                                      startidx2, startidx1, myPeriodicData,
                                                      GFSCAGG!, currfreqG, currk);
                        end
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

            assembleDeltaKI!(DeltaKImat, myDeltaKISystemData, myAllMolData,
                             myPeriodicData, currk);

            
            assembleZIIinvTeeradiation!(totalZIIinvmat, totalTeemat, totalGenvmat,
                                        myAllMolData, myTransData, tt);

            assemblePhimnradiation!(Phimn, myAllMolData, myPeriodicData, totalTeemat,
                                    totalGenvmatwithdiag, totalZIIinvmat, tempmat,
                                    currfreq);

            curroutfilename = string("radiation_", outfilename);
            
            printPhimn(Phimn, curroutfilename, myAllMolData, myTransData, myPeriodicData,
                       tt, currfreq, currk);

            
            assembleZIIinvTeeconduction!(totalZIIinvmat, totalTeemat, DeltaKImat,
                                         myAllMolData, myTransData, tt);

            assemblePhimnconduction!(Phimn, myAllMolData, myPeriodicData, totalTeemat,
                                     DeltaKImat, totalZIIinvmat, tempmat, currfreq);

            curroutfilename = string("conduction_", outfilename);
            
            printPhimn(Phimn, curroutfilename, myAllMolData, myTransData, myPeriodicData,
                       tt, currfreq, currk);

            
            assembleZIIinvTeeboth!(totalZIIinvmat, totalTeemat, totalGenvmat, DeltaKImat,
                                   myAllMolData, myTransData, tt);

            assemblePhimnboth!(Phimn, myAllMolData, myPeriodicData, totalTeemat,
                               totalGenvmat, DeltaKImat, totalZIIinvmat, tempmat, currfreq);

            curroutfilename = string("both_", outfilename);
            
            printPhimn(Phimn, curroutfilename, myAllMolData, myTransData, myPeriodicData,
                       tt, currfreq, currk);
                        

            untransformAtomPos!(myTransData, myAllMolData, tt);

        end
    
    end
  
end

"""

    printPhimn(Phimn::Array{FT, 2}, outfilename::AbstractString,
               myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
               myTransData::TransData{FT},
               myPeriodicData::PeriodicData{FT},
               tt::Integer, currfreq::Union{FT, Complex{FT}},
               currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Print the matrix `Phimn` to the file with name `outfilename`. The elements `Phimn[m, n]`
represent the energy transfer from molecule `m` to molecule `n`, corresponding to the
mathematical quantity ``\\Phi^{(m)}_{n}``. If there is no spatial periodicity,
`Phimn` should be symmetric, and a line of output will be as follows:
`myTransData.transLabels[tt]` `currfreq` `Phimn[1, 1]` `Phimn[1, 2]` ... `Phimn[2, 1]` `Phimn[2, 2]` ...
for example:
`TRANS05 1.2340953873461249e+12 -6.3540963812489034e-01 6.3540963812489034e-01 6.3540963812489034e-01 -6.3540963812489034e-01`  
If there is spatial periodicity, the transpose of `Phimn` at wavevector ``\vec{k}`` is the
equivalent of computing `Phimn` at ``-\vec{k}``, and a line of output will be as follows:
`myTransData.transLabels[tt]` `currfreq` `currk[1]` `currk[2]` `currk[3]` `Phimn[1, 1]` `Phimn[1, 2]` ... `Phimn[2, 1]` `Phimn[2, 2]` ...
for example:
`TRANS05 1.2340953873461249e+12 5.0000000000000000e+07 3.0000000000000000e+07 0.0000000000000000e+00 -6.3540963812489034e-01 6.3540963812489034e-01 6.3540963812489034e-01 -6.3540963812489034e-01`

"""
function printPhimn(Phimn::Array{FT, 2}, outfilename::AbstractString,
                    myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                    myTransData::TransData{FT},
                    myPeriodicData::PeriodicData{FT},
                    tt::Integer, currfreq::Union{FT, Complex{FT}},
                    currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    outfile = open(outfilename, "a+");
    @printf(outfile, "%s %.16e", myTransData.transLabels[tt], abs(currfreq));
    if (myPeriodicData.numdims == 0)
    
        for nn1=1:myAllMolData.nummols
            for nn2=1:myAllMolData.nummols
                @printf(outfile, " %.16e", Phimn[nn1, nn2]);;
            end
        end
    else
        @printf(outfile, " %.16e %.16e %.16e", currk[1], currk[2], currk[3]);
        for nn1=1:myAllMolData.nummols
            for nn2=1:myAllMolData.nummols
                @printf(outfile, " %.16e", Phimn[nn1, nn2]);
            end
        end
        @printf(outfile, "\n%s %.16e %.16e %.16e %.16e",
                myTransData.transLabels[tt], abs(currfreq), -1*currk[1],
                -1*currk[2], -1*currk[3]);
        for nn1=1:myAllMolData.nummols
            for nn2=1:myAllMolData.nummols
                @printf(outfile, " %.16e", Phimn[nn2, nn1]);
            end
        end
    end
    @printf(outfile, "\n");
    close(outfile);
    
end

"""

    assembleZIIinvnoDeltaKI!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                             myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                             myTransData::TransData{FT},
                             tt::Integer)

Assemble the matrix `totalZIIinvmat` in-place for the geometrical transformation with index
`tt`, in the absence of internuclear couplings between different molecules.

"""
function assembleZIIinvnoDeltaKI!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                  myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                  myTransData::TransData{FT},
                                  tt::Integer) where FT<:AbstractFloat

    fill!(totalZIIinvmat, zero(eltype(totalZIIinvmat)));
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

        tempview = view(totalZIIinvmat, startidx:endidx, startidx:endidx);
        tempview .= myAllMolData.uniqueMolDataArray[currumolind].KIBloch;
        addtodiag!(tempview,
                   vecNto3N(1 ./ myAllMolData.uniqueMolDataArray[currumolind].Ye0II));
        tempview .= inv(tempview);
        
        rotate3Nmat!(myAllMolData.uniqueMolDataArray[currumolind].KIBloch,
                     transpose(myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] *
                               myTransData.baselineRot[3*nn-2:3*nn, :]));

    end
  
end

"""

    assembleZIIinvTeeradiation!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                totalTeemat::SharedArray{Complex{FT}, 2},
                                totalGenvmat::SharedArray{Complex{FT}, 2},
                                myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                myTransData::TransData{FT},
                                tt::Integer)

Assemble the matrices `totalZIIinvmat` and `totalTeemat` in-place for the geometrical
transformation with index `tt`, in the absence of internuclear couplings between different
molecules but in the presence of radiative couplings.

"""
function assembleZIIinvTeeradiation!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                     totalTeemat::SharedArray{Complex{FT}, 2},
                                     totalGenvmat::SharedArray{Complex{FT}, 2},
                                     myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                     myTransData::TransData{FT},
                                     tt::Integer) where FT<:AbstractFloat

    fill!(totalZIIinvmat, zero(eltype(totalZIIinvmat)));
    fill!(totalTeemat, zero(eltype(totalTeemat)));

    assembleZIIinvnoDeltaKI!(totalZIIinvmat, myAllMolData, myTransData, tt);
    
    totalTeemat .= -1 .* totalZIIinvmat;

    for nn=1:myAllMolData.nummols

        currumolind = myAllMolData.uniquemolinds[nn];
        if (currumolind < 0)
            currumolind = myAllMolData.dupemolinds[nn];
        end

        startidx = 1;
        if (nn > 1)
            startidx = 3*myAllMolData.cumnumatomslist[nn-1] + 1;
        end
        endidx = 3*myAllMolData.cumnumatomslist[nn];

        tempview = view(totalTeemat, startidx:endidx, startidx:endidx);
        mulMatDiagleft!(tempview,
                        vecNto3N(myAllMolData.uniqueMolDataArray[currumolind].Ke));
        mulMatDiagright!(tempview,
                        vecNto3N(myAllMolData.uniqueMolDataArray[currumolind].Ke));
        addtodiag!(tempview,
                   vecNto3N(1 ./ myAllMolData.uniqueMolDataArray[currumolind].Ye0ee));

    end

    totalTeemat .= inv(totalTeemat .- totalGenvmat);

end

"""

    assembletempmateI!(tempmat::SharedArray{Complex{FT}, 2},
                       totalTeemat::SharedArray{Complex{FT}, 2},
                       totalZIIinvmat::SharedArray{Complex{FT}, 2},
                       myAllMolData::MolSystem{OneMol_WithPhonons{FT}})

Assemble the matrix `tempmat` in-place, equal to
``(Z_{\\mathrm{ee}} - K_{\\mathrm{e}} Z_{\\mathrm{II}}^{-1} K_{\\mathrm{e}})^{-1} K_{\\mathrm{e}} Z_{\\mathrm{II}}^{-1}``.

"""
function assembletempmateI!(tempmat::SharedArray{Complex{FT}, 2},
                            totalTeemat::SharedArray{Complex{FT}, 2},
                            totalZIIinvmat::SharedArray{Complex{FT}, 2},
                            myAllMolData::MolSystem{OneMol_WithPhonons{FT}}) where FT<:AbstractFloat

    tempmat .= totalZIIinvmat;
    for nn=1:myAllMolData.nummols

        currumolind = myAllMolData.uniquemolinds[nn];
        if (currumolind < 0)
            currumolind = myAllMolData.dupemolinds[nn];
        end

        startidx = 1;
        if (nn > 1)
            startidx = 3*myAllMolData.cumnumatomslist[nn-1] + 1;
        end
        endidx = 3*myAllMolData.cumnumatomslist[nn];

        tempview = view(tempmat, startidx:endidx, 1:3*myAllMolData.cumnumatomslist[end]);
        mulMatDiagleft!(tempview,
                        vecNto3N(myAllMolData.uniqueMolDataArray[currumolind].Ke));

    end
    tempmat .= totalTeemat * tempmat; 

end

"""

    assemblePhimnwithBe!(Phimn::Array{FT, 2},
                         myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                         myPeriodicData::PeriodicData{FT},
                         tempmat::SharedArray{Complex{FT}, 2},
                         DeltaZImat::SharedArray{Complex{FT}, 2},
                         currfreq::Union{FT, Complex{FT}})

Assemble the contributions to the matrix `Phimn` in-place, corresponding to the fluctuating
electronic sources (accounting for absorption in both electrons due to radiation and nuclei
due to conduction); contributions are added to the existing values of `Phimn`.

"""
function assemblePhimnwithBe!(Phimn::Array{FT, 2},
                              myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                              myPeriodicData::PeriodicData{FT},
                              tempmat::SharedArray{Complex{FT}, 2},
                              DeltaZImat::SharedArray{Complex{FT}, 2},
                              currfreq::Union{FT, Complex{FT}}) where FT<:AbstractFloat

    currPhi = zero(eltype(Phimn));
  
    for nn2=1:myAllMolData.nummols

        startidx2 = 1;
        if (nn2 > 1)
            startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
        end
        endidx2 = 3*myAllMolData.cumnumatomslist[nn2];
        DeltaZIview = view(DeltaZImat, startidx2:endidx2,
                           1:3*myAllMolData.cumnumatomslist[end]);

        startmol1 = 1;
        if (myPeriodicData.numdims == 0)
            startmol1 = nn2;
        end
        
        for nn1=startmol1:myAllMolData.nummols
          
            startidx1 = 1;
            if (nn1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
            end
            endidx1 = 3*myAllMolData.cumnumatomslist[nn1];
            
            Tview = view(tempmat, 1:3*myAllMolData.cumnumatomslist[end], startidx1:endidx1);
            
            currmolind1 = myAllMolData.uniquemolinds[nn1];
            if (currmolind1 < 0)
                currmolind1 = myAllMolData.dupemolinds[nn1];
            end

            currPhi = currfreq *
            imag(tr(mulMatDiagleft(adjoint(Tview)[:, startidx2:endidx2] * DeltaZIview *
                                   Tview,
                                   vecNto3N(myAllMolData.uniqueMolDataArray[currmolind1].Be))));

            Phimn[nn1, nn2] += currPhi;
            
            if (myPeriodicData.numdims == 0 && nn1 != nn2)
                Phimn[nn2, nn1] += currPhi;
            end
        end
        
    end

end

"""

    assemblePhimnwithBI!(Phimn::Array{FT, 2},
                         myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                         myPeriodicData::PeriodicData{FT},
                         tempmat::SharedArray{Complex{FT}, 2},
                         DeltaZImat::SharedArray{Complex{FT}, 2},
                         currfreq::Union{FT, Complex{FT}})

Assemble the contributions to the matrix `Phimn` in-place, corresponding to the fluctuating
nuclear sources (accounting for absorption in both electrons due to radiation and nuclei
due to conduction); contributions are added to the existing values of `Phimn`.

"""
function assemblePhimnwithBI!(Phimn::Array{FT, 2},
                              myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                              myPeriodicData::PeriodicData{FT},
                              tempmat::SharedArray{Complex{FT}, 2},
                              DeltaZImat::SharedArray{Complex{FT}, 2},
                              currfreq::Union{FT, Complex{FT}}) where FT<:AbstractFloat

    currPhi = zero(eltype(Phimn));

    for nn2=1:myAllMolData.nummols

        startidx2 = 1;
        if (nn2 > 1)
            startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
        end
        endidx2 = 3*myAllMolData.cumnumatomslist[nn2];
        DeltaZIview = view(DeltaZImat, startidx2:endidx2,
                           1:3*myAllMolData.cumnumatomslist[end]);

        startmol1 = 1;
        if (myPeriodicData.numdims == 0)
            startmol1 = nn2;
        end
        
        for nn1=startmol1:myAllMolData.nummols
          
            startidx1 = 1;
            if (nn1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
            end
            endidx1 = 3*myAllMolData.cumnumatomslist[nn1];
            
            Tview = view(tempmat, 1:3*myAllMolData.cumnumatomslist[end], startidx1:endidx1);
            
            currmolind1 = myAllMolData.uniquemolinds[nn1];
            if (currmolind1 < 0)
                currmolind1 = myAllMolData.dupemolinds[nn1];
            end

            currPhi = currfreq *
            imag(tr(mulMatDiagleft(adjoint(Tview)[:, startidx2:endidx2] * DeltaZIview *
                                   Tview,
                                   vecNto3N(myAllMolData.uniqueMolDataArray[currmolind1].BI))));

            Phimn[nn1, nn2] += currPhi;
            
            if (myPeriodicData.numdims == 0 && nn1 != nn2)
                Phimn[nn2, nn1] += currPhi;
            end
        end
        
    end

end

"""

    assemblePhimnradiation!(Phimn::Array{FT, 2},
                            myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                            myPeriodicData::PeriodicData{FT},
                            totalTeemat::SharedArray{Complex{FT}, 2},
                            totalGenvmatwithdiag::SharedArray{Complex{FT}, 2},
                            totalZIIinvmat::SharedArray{Complex{FT}, 2},
                            tempmat::SharedArray{Complex{FT}, 2},
                            currfreq::Union{FT, Complex{FT}})

Assemble the matrix `Phimn` in-place for purely radiative couplings.

"""
function assemblePhimnradiation!(Phimn::Array{FT, 2},
                                 myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                 myPeriodicData::PeriodicData{FT},
                                 totalTeemat::SharedArray{Complex{FT}, 2},
                                 totalGenvmatwithdiag::SharedArray{Complex{FT}, 2},
                                 totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                 tempmat::SharedArray{Complex{FT}, 2},
                                 currfreq::Union{FT, Complex{FT}}) where FT<:AbstractFloat

    fill!(tempmat, zero(eltype(tempmat)));
  
    fill!(Phimn, zero(eltype(Phimn)));

    ## first: radiation into molecule n due to electronic sources in molecule m 
    tempmat .= totalTeemat;

    lmul!(-1, totalGenvmatwithdiag);
    assemblePhimnwithBe!(Phimn, myAllMolData, myPeriodicData, tempmat, totalGenvmatwithdiag,
                         currfreq);
    lmul!(-1, totalGenvmatwithdiag);

    ## second: radiation into molecule n due to nuclear sources in molecule m 
    assembletempmateI!(tempmat, totalTeemat, totalZIIinvmat, myAllMolData);

    lmul!(-1, totalGenvmatwithdiag);
    assemblePhimnwithBI!(Phimn, myAllMolData, myPeriodicData, tempmat, totalGenvmatwithdiag,
                         currfreq);
    lmul!(-1, totalGenvmatwithdiag);
    
end

"""

    assembleZIIinvwithDeltaKI!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                               DeltaKImat::SharedArray{Complex{FT}, 2},
                               myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                               myTransData::TransData{FT},
                               tt::Integer)

Assemble the matrix `totalZIIinvmat` in-place for the geometrical transformation with index
`tt`, in the presence of internuclear couplings between different molecules.

"""
function assembleZIIinvwithDeltaKI!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                    DeltaKImat::SharedArray{Complex{FT}, 2},
                                    myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                    myTransData::TransData{FT},
                                    tt::Integer) where FT<:AbstractFloat

    fill!(totalZIIinvmat, zero(eltype(totalZIIinvmat)));
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

        tempview = view(totalZIIinvmat, startidx:endidx, startidx:endidx);
        tempview .= myAllMolData.uniqueMolDataArray[currumolind].KIBloch;
        addtodiag!(tempview,
                   vecNto3N(1 ./ myAllMolData.uniqueMolDataArray[currumolind].Ye0II));

        rotate3Nmat!(myAllMolData.uniqueMolDataArray[currumolind].KIBloch,
                     transpose(myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] *
                               myTransData.baselineRot[3*nn-2:3*nn, :]));

    end

    totalZIIinvmat .= inv(totalZIIinvmat .+ DeltaKImat);
  
end

"""

    assembleZIIinvTeeconduction!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                 totalTeemat::SharedArray{Complex{FT}, 2},
                                 DeltaKImat::SharedArray{Complex{FT}, 2},
                                 myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                 myTransData::TransData{FT},
                                 tt::Integer)

Assemble the matrices `totalZIIinvmat` and `totalTeemat` in-place for the geometrical
transformation with index `tt`, in the presence of internuclear couplings between different
molecules but in the absence of radiative couplings.

"""
function assembleZIIinvTeeconduction!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                      totalTeemat::SharedArray{Complex{FT}, 2},
                                      DeltaKImat::SharedArray{Complex{FT}, 2},
                                      myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                      myTransData::TransData{FT},
                                      tt::Integer) where FT<:AbstractFloat

    fill!(totalZIIinvmat, zero(eltype(totalZIIinvmat)));
    fill!(totalTeemat, zero(eltype(totalTeemat)));

    assembleZIIinvwithDeltaKI!(totalZIIinvmat, DeltaKImat, myAllMolData, myTransData, tt);

    totalTeemat .= -1 .* totalZIIinvmat;

    for nn2=1:myAllMolData.nummols

        currmolind2 = myAllMolData.uniquemolinds[nn2];
        if (currmolind2 < 0)
            currmolind2 = myAllMolData.dupemolinds[nn2];
        end

        startidx2 = 1;
        if (nn2 > 1)
            startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
        end
        endidx2 = 3*myAllMolData.cumnumatomslist[nn2];

        for nn1=1:myAllMolData.nummols

            currmolind1 = myAllMolData.uniquemolinds[nn1];
            if (currmolind1 < 0)
                currmolind1 = myAllMolData.dupemolinds[nn1];
            end

            startidx1 = 1;
            if (nn1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
            end
            endidx1 = 3*myAllMolData.cumnumatomslist[nn1];

            tempview = view(totalTeemat, startidx1:endidx1, startidx2:endidx2);
            mulMatDiagleft!(tempview,
                            vecNto3N(myAllMolData.uniqueMolDataArray[currmolind1].Ke));
            mulMatDiagright!(tempview,
                             vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].Ke));

        end

        tempview = view(totalTeemat, startidx2:endidx2, startidx2:endidx2)
        addtodiag!(tempview,
                   vecNto3N(1 ./ myAllMolData.uniqueMolDataArray[currmolind2].Ye0ee));

    end

    totalTeemat .= inv(totalTeemat);

end

"""

    assembletempmatIe!(tempmat::SharedArray{Complex{FT}, 2},
                       totalTeemat::SharedArray{Complex{FT}, 2},
                       totalZIIinvmat::SharedArray{Complex{FT}, 2},
                       myAllMolData::MolSystem{OneMol_WithPhonons{FT}})

Assemble the matrix `tempmat` in-place, equal to
``Z_{\\mathrm{II}}^{-1} K_{\\mathrm{e}} (Z_{\\mathrm{ee}} - K_{\\mathrm{e}} Z_{\\mathrm{II}}^{-1} K_{\\mathrm{e}})^{-1}``.

"""
function assembletempmatIe!(tempmat::SharedArray{Complex{FT}, 2},
                            totalTeemat::SharedArray{Complex{FT}, 2},
                            totalZIIinvmat::SharedArray{Complex{FT}, 2},
                            myAllMolData::MolSystem{OneMol_WithPhonons{FT}}) where FT<:AbstractFloat

    tempmat .= totalZIIinvmat;
    for nn=1:myAllMolData.nummols

        currumolind = myAllMolData.uniquemolinds[nn];
        if (currumolind < 0)
            currumolind = myAllMolData.dupemolinds[nn];
        end

        startidx = 1;
        if (nn > 1)
            startidx = 3*myAllMolData.cumnumatomslist[nn-1] + 1;
        end
        endidx = 3*myAllMolData.cumnumatomslist[nn];

        tempview = view(tempmat, 1:3*myAllMolData.cumnumatomslist[end], startidx:endidx);
        mulMatDiagright!(tempview,
                         vecNto3N(myAllMolData.uniqueMolDataArray[currumolind].Ke));

    end
    tempmat .= totalTeemat * tempmat; 

end

"""

    assembletempmatII!(tempmat::SharedArray{Complex{FT}, 2},
                       totalTeemat::SharedArray{Complex{FT}, 2},
                       totalZIIinvmat::SharedArray{Complex{FT}, 2},
                       myAllMolData::MolSystem{OneMol_WithPhonons{FT}})

Assemble the matrix `tempmat` in-place, equal to
``Z_{\\mathrm{II}}^{-1} + Z_{\\mathrm{II}}^{-1} K_{\\mathrm{e}} (Z_{\\mathrm{ee}} - K_{\\mathrm{e}} Z_{\\mathrm{II}}^{-1} K_{\\mathrm{e}})^{-1} K_{\\mathrm{e}} Z_{\\mathrm{II}}^{-1}``.

"""
function assembletempmatII!(tempmat::SharedArray{Complex{FT}, 2},
                            totalTeemat::SharedArray{Complex{FT}, 2},
                            totalZIIinvmat::SharedArray{Complex{FT}, 2},
                            myAllMolData::MolSystem{OneMol_WithPhonons{FT}}) where FT<:AbstractFloat

    tempmat .= totalTeemat;
    for nn2=1:myAllMolData.nummols

        currmolind2 = myAllMolData.uniquemolinds[nn2];
        if (currmolind2 < 0)
            currmolind2 = myAllMolData.dupemolinds[nn2];
        end

        startidx2 = 1;
        if (nn2 > 1)
            startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
        end
        endidx2 = 3*myAllMolData.cumnumatomslist[nn2];

        for nn1=1:myAllMolData.nummols

            currmolind1 = myAllMolData.uniquemolinds[nn1];
            if (currmolind1 < 0)
                currmolind1 = myAllMolData.dupemolinds[nn1];
            end

            startidx1 = 1;
            if (nn1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
            end
            endidx1 = 3*myAllMolData.cumnumatomslist[nn1];

            tempview = view(tempmat, startidx1:endidx1, startidx2:endidx2);
            mulMatDiagleft!(tempview,
                            vecNto3N(myAllMolData.uniqueMolDataArray[currmolind1].Ke));
            mulMatDiagright!(tempview,
                             vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].Ke));

        end

    end
    tempmat .= totalZIIinvmat .+ (totalZIIinvmat * tempmat * totalZIIinvmat); 

end

"""

    assemblePhimnconduction!(Phimn::Array{FT, 2},
                             myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                             myPeriodicData::PeriodicData{FT},
                             totalTeemat::SharedArray{Complex{FT}, 2},
                             DeltaKImat::SharedArray{Complex{FT}, 2},
                             totalZIIinvmat::SharedArray{Complex{FT}, 2},
                             tempmat::SharedArray{Complex{FT}, 2},
                             currfreq::Union{FT, Complex{FT}})

Assemble the matrix `Phimn` in-place for purely conductive couplings.

"""
function assemblePhimnconduction!(Phimn::Array{FT, 2},
                                  myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                  myPeriodicData::PeriodicData{FT},
                                  totalTeemat::SharedArray{Complex{FT}, 2},
                                  DeltaKImat::SharedArray{Complex{FT}, 2},
                                  totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                  tempmat::SharedArray{Complex{FT}, 2},
                                  currfreq::Union{FT, Complex{FT}}) where FT<:AbstractFloat

    fill!(tempmat, zero(eltype(tempmat)));

    fill!(Phimn, zero(eltype(Phimn)));

    ## first: conduction into molecule n due to electronic sources in molecule m 
    assembletempmatIe!(tempmat, totalTeemat, totalZIIinvmat, myAllMolData);

    assemblePhimnwithBe!(Phimn, myAllMolData, myPeriodicData, tempmat, DeltaKImat,
                         currfreq);

    ## second: conduction into molecule n due to nuclear sources in molecule m 
    assembletempmatII!(tempmat, totalTeemat, totalZIIinvmat, myAllMolData);    

    assemblePhimnwithBI!(Phimn, myAllMolData, myPeriodicData, tempmat, DeltaKImat,
                         currfreq);
    
end

"""

    assembleZIIinvTeeboth!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                           totalTeemat::SharedArray{Complex{FT}, 2},
                           totalGenvmat::SharedArray{Complex{FT}, 2},
                           DeltaKImat::SharedArray{Complex{FT}, 2},
                           myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                           myTransData::TransData{FT},
                           tt::Integer)

Assemble the matrices `totalZIIinvmat` and `totalTeemat` in-place for the geometrical
transformation with index `tt`, in the presence of both radiative and internuclear couplings
between different molecules.

"""
function assembleZIIinvTeeboth!(totalZIIinvmat::SharedArray{Complex{FT}, 2},
                                totalTeemat::SharedArray{Complex{FT}, 2},
                                totalGenvmat::SharedArray{Complex{FT}, 2},
                                DeltaKImat::SharedArray{Complex{FT}, 2},
                                myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                myTransData::TransData{FT},
                                tt::Integer) where FT<:AbstractFloat

    fill!(totalZIIinvmat, zero(eltype(totalZIIinvmat)));
    fill!(totalTeemat, zero(eltype(totalTeemat)));

    assembleZIIinvwithDeltaKI!(totalZIIinvmat, DeltaKImat, myAllMolData, myTransData, tt);

    totalTeemat .= -1 .* totalZIIinvmat;

    for nn2=1:myAllMolData.nummols

        currmolind2 = myAllMolData.uniquemolinds[nn2];
        if (currmolind2 < 0)
            currmolind2 = myAllMolData.dupemolinds[nn2];
        end

        startidx2 = 1;
        if (nn2 > 1)
            startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
        end
        endidx2 = 3*myAllMolData.cumnumatomslist[nn2];

        for nn1=1:myAllMolData.nummols

            currmolind1 = myAllMolData.uniquemolinds[nn1];
            if (currmolind1 < 0)
                currmolind1 = myAllMolData.dupemolinds[nn1];
            end

            startidx1 = 1;
            if (nn1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
            end
            endidx1 = 3*myAllMolData.cumnumatomslist[nn1];

            tempview = view(totalTeemat, startidx1:endidx1, startidx2:endidx2);
            mulMatDiagleft!(tempview,
                            vecNto3N(myAllMolData.uniqueMolDataArray[currmolind1].Ke));
            mulMatDiagright!(tempview,
                             vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].Ke));

        end

        tempview = view(totalTeemat, startidx2:endidx2, startidx2:endidx2)
        addtodiag!(tempview,
                   vecNto3N(1 ./ myAllMolData.uniqueMolDataArray[currmolind2].Ye0ee));

    end

    totalTeemat .= inv(totalTeemat .- totalGenvmat);

end

"""

    assemblePhimnboth!(Phimn::Array{FT, 2},
                       myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                       myPeriodicData::PeriodicData{FT},
                       totalTeemat::SharedArray{Complex{FT}, 2},
                       totalGenvmatwithdiag::SharedArray{Complex{FT}, 2},
                       DeltaKImat::SharedArray{Complex{FT}, 2},
                       totalZIIinvmat::SharedArray{Complex{FT}, 2},
                       tempmat::SharedArray{Complex{FT}, 2},
                       currfreq::Union{FT, Complex{FT}})

Assemble the matrix `Phimn` in-place for both radiative and conductive couplings.

"""
function assemblePhimnboth!(Phimn::Array{FT, 2},
                            myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                            myPeriodicData::PeriodicData{FT},
                            totalTeemat::SharedArray{Complex{FT}, 2},
                            totalGenvmatwithdiag::SharedArray{Complex{FT}, 2},
                            DeltaKImat::SharedArray{Complex{FT}, 2},
                            totalZIIinvmat::SharedArray{Complex{FT}, 2},
                            tempmat::SharedArray{Complex{FT}, 2},
                            currfreq::Union{FT, Complex{FT}}) where FT<:AbstractFloat

    fill!(tempmat, zero(eltype(tempmat)));

    fill!(Phimn, zero(eltype(Phimn)));

    ## first: radiation into molecule n due to electronic sources in molecule m 
    tempmat .= totalTeemat;

    lmul!(-1, totalGenvmatwithdiag);
    assemblePhimnwithBe!(Phimn, myAllMolData, myPeriodicData, tempmat, totalGenvmatwithdiag,
                         currfreq);
    lmul!(-1, totalGenvmatwithdiag);

    ## second: radiation into molecule n due to nuclear sources in molecule m 
    assembletempmateI!(tempmat, totalTeemat, totalZIIinvmat, myAllMolData);

    lmul!(-1, totalGenvmatwithdiag);
    assemblePhimnwithBI!(Phimn, myAllMolData, myPeriodicData, tempmat, totalGenvmatwithdiag,
                         currfreq);
    lmul!(-1, totalGenvmatwithdiag);
    
    ## third: conduction into molecule n due to electronic sources in molecule m 
    assembletempmatIe!(tempmat, totalTeemat, totalZIIinvmat, myAllMolData);

    assemblePhimnwithBe!(Phimn, myAllMolData, myPeriodicData, tempmat, DeltaKImat,
                         currfreq);

    ## fourth: conduction into molecule n due to nuclear sources in molecule m 
    assembletempmatII!(tempmat, totalTeemat, totalZIIinvmat, myAllMolData);    

    assemblePhimnwithBI!(Phimn, myAllMolData, myPeriodicData, tempmat, DeltaKImat,
                         currfreq);
    
end

"""

    assembleGenvDiagBlock!(totalGenvmat::SharedArray{Complex{FT}, 2},
                           totalGenvmatwithdiag::SharedArray{Complex{FT}, 2},
                           myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                           myTransData::TransData{FT}, nn2::Integer,
                           currmolind2::Integer, startidx2::Integer,
                           tt::Integer, myPeriodicData::PeriodicData{FT},
                           GFSCAGG!::Function,
                           currfreqG::Union{FT, Complex{FT}},
                           currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the molecular diagonal block of the scattering Green's function interaction matrix
of molecule `nn2` (whose unique molecular index `currmolind2` is precomputed and passed as
an argument) at frequency `currfreqG` and wavevector `currk` using `GFSCAGG!`, then add the
vacuum contribution of molecule `nn2`, performing any required matrix rotations
corresponding to transformation `tt` for molecule `nn2`. Stamp that into
`totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`,
and added to coincident contributions into
`totalGenvmatwithdiag[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
Finally, undo any rotation corresponding to transformation `tt` in molecule `nn2` for reuse
with future transformations. This assumes the presence of phonons.

"""
function assembleGenvDiagBlock!(totalGenvmat::SharedArray{Complex{FT}, 2},
                                totalGenvmatwithdiag::SharedArray{Complex{FT}, 2},
                                myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                myTransData::TransData{FT}, nn2::Integer,
                                currmolind2::Integer, startidx2::Integer,
                                tt::Integer, myPeriodicData::PeriodicData{FT},
                                GFSCAGG!::Function,
                                currfreqG::Union{FT, Complex{FT}},
                                currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalGenvmatwithdiag, (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));

    fill!(tempview, zero(eltype(tempview)));

    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn2, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                 myTransData.baselineRot[3*nn2-2:3*nn2, :]);
    
    tempview .+= myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf[:, :];

    totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                 (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])] .=
    tempview;

    addtodiag!(tempview,
               vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinfdiag));
    
    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 transpose(myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                           myTransData.baselineRot[3*nn2-2:3*nn2, :]));

end

"""

    assembleGenvOffDiagBlock!(totalGenvmat::SharedArray{Complex{FT}, 2},
                              totalGenvmatwithdiag::SharedArray{Complex{FT}, 2},
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
the left by the electronic polarizabilities, add the vacuum contribution, and stamp that
into
`totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`
and
`totalGenvmatwithdiag[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
If there is no spatial periodicity, do the same for the (unconjugated) transpose, swapping
`nn2` and `nn1` without performing further Green's function computations.

"""
function assembleGenvOffDiagBlock!(totalGenvmat::SharedArray{Complex{FT}, 2},
                                   totalGenvmatwithdiag::SharedArray{Complex{FT}, 2},
                                   myAllMolData::MolSystem{<:OneMol{FT}},
                                   myTransData::TransData{FT},
                                   nn1::Integer, nn2::Integer,
                                   startidx1::Integer,
                                   startidx2::Integer,
                                   myPeriodicData::PeriodicData{FT},
                                   GFSCAGG!::Function,
                                   currfreqG::Union{FT, Complex{FT}},
                                   currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalGenvmatwithdiag,
                    (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));
    fill!(tempview, zero(eltype(tempview)));

    assembleGvacMolBlock!(tempview, myAllMolData, myTransData, nn1, nn2, 1, 1,
                          myPeriodicData, currfreqG, currk);
    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn1, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);

    totalGenvmat[(startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]),
                         (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])] .=
    tempview;

    
    if (myPeriodicData.numdims == 0)
        totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                     (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1])] .=
        transpose(deepcopy(tempview));
        totalGenvmatwithdiag[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                             (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1])] .=
        totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                     (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1])]; 
   end

end
