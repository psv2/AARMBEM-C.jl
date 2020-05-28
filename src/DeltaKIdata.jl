"""

    DeltaKIMolBlockData{FT<:AbstractFloat}

Represent the internuclear couplings between a pair of molecules, assumed to be computed via
a Morse potential, with the possibility of Bloch periodicity.

See also: [`DeltaKISystemData{FT}`](@ref)

"""
struct DeltaKIMolBlockData{FT<:AbstractFloat}

    mol1::Integer
    mol2::Integer
    k0molblock::Array{FT, 2}
    a0molblock::Array{FT, 2}
    Einfmolblock::Array{FT, 2}
    DeltaKImolblock::SharedArray{FT, 2}
    
    numblocks1::Integer
    blocklist1::UnitRange{<:Integer}
    numblocks2::Integer
    blocklist2::UnitRange{<:Integer}
    
end

"""

    DeltaKISystemData{FT<:AbstractFloat}

Represent the internuclear couplings between every pair of molecules, assumed to be computed
via a Morse potential, with the possibility of Bloch periodicity, and represent the diagonal
elements of that internuclear spring constant matrix for the full system based on the
off-diagonal blocks.

See also: [`DeltaKIMolBlockData{FT}`](@ref)

"""
struct DeltaKISystemData{FT<:AbstractFloat}

    DeltaKImolblockdataarray::Array{DeltaKIMolBlockData{FT}, 1}
    DeltaKIdiagelements::Array{FT, 1}
    numcouplings::Integer

end

"""

    keffMorse(r::FT, k0::FT, a0::FT, Einf::FT) where FT<:AbstractFloat

Compute and return the effective Morse spring constant between two particles at distance `r`
(assumed to be nonnegative), with equilibrium spring constant `k0` and separation `a0` and
with dissociation energy `Einf`. If ``U(r)`` is the potential function, then this returns
the quantity ``(1/(r - a_{0}))\\partial U/\\partial r``.

"""
function keffMorse(r::FT, k0::FT, a0::FT, Einf::FT) where FT<:AbstractFloat

    if (k0 == zero(k0) || a0 == zero(a0) || Einf == zero(Einf))
        return zero(r);
    elseif (r == a0)
        return -1*k0;
    end

    x = sqrt(k0/(2*Einf))*(r - a0);
    return k0*exp(-1*x)*expm1(-1*x)/x;

end

"""

    readDeltaKISystemData(mollistfilename::AbstractString,
                          myAllMolData::MolSystem{OneMol_WithPhonons{FT}}) where FT<:AbstractFloat

Read a configuration file `mollistfilename` (a string) to extract data detailing which pairs
of molecules have nontrivial internuclear couplings between them, handling duplicate entries
appropriately.

"""
function readDeltaKISystemData(mollistfilename::AbstractString,
                               myAllMolData::MolSystem{OneMol_WithPhonons{FT}}) where FT<:AbstractFloat

    molstrarray = strip.(chomp.(readlines(mollistfilename)));
    couplingstartinds = findall(map(u -> startswith(u, "COUPLING"), molstrarray));
    couplingendinds = findall(molstrarray .== "ENDCOUPLING");
    numcouplings = length(couplingendinds);
    if (length(couplingstartinds) != numcouplings)
        error(string(mollistfilename, ": number of COUPLING and ENDCOUPLING statements",
                     " should be the same (may be zero)"));
    end

    # if there are no couplings, return a structure with an empty array
    if (numcouplings == 0)
        return DeltaKISystemData{FT}(Array{DeltaKIMolBlockData{FT}, 1}(undef, 0),
                                     zeros(FT, myAllMolData.cumnumatomslist[end]),
                                     numcouplings);
    end

    DeltaKImolblockdataarray = Array{DeltaKIMolBlockData{FT}, 1}(undef, numcouplings);
    molpairarray = Array{Integer, 2}(undef, numcouplings, 2);
    for nn=1:numcouplings

        currcouplingstrarray = strip.(chomp.(molstrarray[couplingstartinds[nn]+1:couplingendinds[nn]-1]));
        DeltaKImolblockdataarray[nn] = readDeltaKIBlockData(currcouplingstrarray,
                                                            myAllMolData, mollistfilename);
        molpairarray[nn, 1] = DeltaKImolblockdataarray[nn].mol1;
        molpairarray[nn, 2] = DeltaKImolblockdataarray[nn].mol2;

    end
    
    # need to figure out which blocks have been duplicated or have self-couplings,
    # in order to remove them
    rowisdupe = falses(numcouplings);

    for nn1=2:numcouplings

        for nn2=1:(nn1-1)

            if (molpairarray[nn1, 1] == molpairarray[nn2, 1] &&
                molpairarray[nn1, 2] == molpairarray[nn2, 2])

                rowisdupe[nn1] = true;
                break;

            end

        end

    end
    deleteidxarr = findall(rowisdupe .| (molpairarray[:, 1] .== molpairarray[:, 2]));
    if (length(deleteidxarr) > 0)

        deleteat!(DeltaKImolblockdataarray, deleteidxarr);
        println(string("If ", mollistfilename, " has multiple consecutive pairs of ",
                       "COUPLING and ENDCOUPLING statements, and there are multiple such ",
                       "pairs where the lines of the form: \"MOL1 \${mol1label}\" and ",
                       "\"MOL2 \${mol2label}\", have the same (or transposed) molecular ",
                       "label pairs mol1label and mol2label, only the first occurrence ",
                       "will be used, even if later occurrences specify different input ",
                       "data files; additionally, any self-coupling blocks specifying ",
                       "\${mol1label} to be the same as \${mol2label} will be ignored"));

    end

    return DeltaKISystemData{FT}(DeltaKImolblockdataarray,
                                 zeros(FT, myAllMolData.cumnumatomslist[end]),
                                 length(DeltaKImolblockdataarray));

end

function readDeltaKIBlockData(couplingstrarray::Array{<:AbstractString, 1},
                              myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                              mollistfilename::AbstractString) where FT<:AbstractFloat

    mol1, mol2, flipmol = extractMolPair(couplingstrarray, myAllMolData, mollistfilename);
    k0, a0, Einf = extractGeneralCouplingData(couplingstrarray, mollistfilename, FT)
    if (flipmol)
        k0 = transpose(k0) .+ zero(eltype(k0));
        a0 = transpose(a0) .+ zero(eltype(a0));
        Einf = transpose(Einf) .+ zero(eltype(Einf));
    end

    numblocks1 = div(size(k0, 1), myAllMolData.numatomslist[mol1]);
    numblocks2 = div(size(k0, 2), myAllMolData.numatomslist[mol2]);
    blocklist1 = div(1 - numblocks1, 2):div(numblocks1 - 1, 2);
    blocklist2 = div(1 - numblocks2, 2):div(numblocks2 - 1, 2);

    return DeltaKIMolBlockData{FT}(mol1, mol2, k0, a0, Einf, zero(k0),
                                   numblocks1, blocklist1, numblocks2, blocklist2);

end

function extractMolPair(couplingstrarray::Array{<:AbstractString, 1},
                        myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                        mollistfilename::AbstractString) where FT<:AbstractFloat

    curridx1 = findfirst(map(u -> startswith(u, "MOL1 "), couplingstrarray));
    if (curridx1 == nothing)
        error(string("If ", mollistfilename, " has at least one consecutive pair of ",
                     "COUPLING and ENDCOUPLING statements, then it must contain a line of ",
                     "the form: MOL1 \${mol1label}"));
    end
    curridx2 = findfirst(map(u -> startswith(u, "MOL2 "), couplingstrarray));
    if (curridx2 == nothing)
        error(string("If ", mollistfilename, " has at least one consecutive pair of ",
                     "COUPLING and ENDCOUPLING statements, then it must contain a line of ",
                     "the form: MOL2 \${mol2label}"));
    end

    mol1label = replace(couplingstrarray[curridx1], "MOL1 " => "");
    mol2label = replace(couplingstrarray[curridx2], "MOL2 " => "");

    mol1 = findfirst(myAllMolData.molLabels .== mol1label);
    if (mol1 == nothing)
        error(string("If ", mollistfilename, " has at least one consecutive pair of ",
                     "COUPLING and ENDCOUPLING statements, then in the line of ",
                     "the form: \"MOL1 \${mol1label}\", mol1label must correspond to a ",
                     "valid molecular label from elsewhere in the same file"));
    end
    mol2 = findfirst(myAllMolData.molLabels .== mol2label);
    if (mol2 == nothing)
        error(string("If ", mollistfilename, " has at least one consecutive pair of ",
                     "COUPLING and ENDCOUPLING statements, then in the line of ",
                     "the form: \"MOL2 \${mol2label}\", mol1label must correspond to a ",
                     "valid molecular label from elsewhere in the same file"));
    end

    flipmol = (mol2 < mol1);
    if (mol2 == mol1)
       error(string("If ", mollistfilename, " has at least one consecutive pair of ",
                     "COUPLING and ENDCOUPLING statements, then in the lines of ",
                     "the form: \"MOL1 \${mol1label}\" and \"MOL2 \${mol2label}\", ",
                     "mol1label and mol2label must be distinct molecular labels"));
    end
    if (flipmol)
        println(string("If ", mollistfilename, " has at least one consecutive pair of ",
                     "COUPLING and ENDCOUPLING statements, then in the lines of ",
                     "the form: \"MOL1 \${mol1label}\" and \"MOL2 \${mol2label}\", ",
                     "it is customary that mol1label refer to a molecular label that ",
                     "comes before mol2label in the same file; all data matrices will ",
                     "therefore be transposed"));
    end

    return mol1, mol2, flipmol;

end

function extractGeneralCouplingData(couplingstrarray::Array{<:AbstractString, 1},
                                    mollistfilename::AbstractString, T::DataType=Float64)

    k0 = extractCouplingDataMat("k0", couplingstrarray, mollistfilename, T);
    a0 = extractCouplingDataMat("a0", couplingstrarray, mollistfilename, T);
    Einf = extractCouplingDataMat("Einf", couplingstrarray, mollistfilename, T);
    return k0, a0, Einf;

end

function extractCouplingDataMat(datalabel::AbstractString,
                                couplingstrarray::Array{<:AbstractString, 1},
                                mollistfilename::AbstractString, T::DataType=Float64)

    curridx = findfirst(map(u -> startswith(u, string(datalabel, " ")), couplingstrarray));
    if (curridx == nothing)
        error(string("If ", mollistfilename, " has at least one consecutive pair of ",
                     "COUPLING and ENDCOUPLING statements, then it must contain a line of ",
                     "the form: ", datalabel, "\${", datalabel, "filename}"));
    end
    return readdlm(replace(couplingstrarray[curridx], string(datalabel, " ") => ""), ',', T);

end
    
function assembleDeltaKImolblock!(myDeltaKImolblockdata::DeltaKIMolBlockData{FT},
                                  myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                  myTransData::TransData{FT},
                                  myPeriodicData::PeriodicData{FT}) where FT<:AbstractFloat

    mymol1 = myDeltaKImolblockdata.mol1;
    mymol2 = myDeltaKImolblockdata.mol2;

    mynumblocks1 = myDeltaKImolblockdata.numblocks1;
    mynumblocks2 = myDeltaKImolblockdata.numblocks2;
    myblocklist1 = myDeltaKImolblockdata.blocklist1;
    myblocklist2 = myDeltaKImolblockdata.blocklist2;
    
    mynumatoms1 = myAllMolData.numatomslist[mymol1];
    mynumatoms2 = myAllMolData.numatomslist[mymol2];
    
    fill!(myDeltaKImolblockdata.DeltaKImolblock, zero(FT));
    atomposqq = zeros(FT, 3);
    atompospp = zeros(FT, 3);
    mylatticevec = zeros(FT, 3);
    
    for nn2=1:mynumblocks2
        for nn1=1:mynumblocks1

            mylatticevec .= myblocklist1[nn1] .* myPeriodicData.latticevecs[:,1] .+
            myblocklist2[nn2] .* myPeriodicData.latticevecs[nn2];

            A = view(myDeltaKImolblockdata.DeltaKImolblock,
                     (nn1-1)*mynumatoms1+1:nn1*mynumatoms1,
                     (nn2-1)*mynumatoms2+1:nn2*mynumatoms2);

            k0 = view(myDeltaKImolblockdata.k0molblock,
                      (nn1-1)*mynumatoms1+1:nn1*mynumatoms1,
                      (nn2-1)*mynumatoms2+1:nn2*mynumatoms2);
            
            a0 = view(myDeltaKImolblockdata.a0molblock,
                      (nn1-1)*mynumatoms1+1:nn1*mynumatoms1,
                      (nn2-1)*mynumatoms2+1:nn2*mynumatoms2);
            
            Einf = view(myDeltaKImolblockdata.Einfmolblock,
                        (nn1-1)*mynumatoms1+1:nn1*mynumatoms1,
                        (nn2-1)*mynumatoms2+1:nn2*mynumatoms2);

            for qq=1:mynumatoms2

                atomposqq .= myTransData.allAtomPos[:, qq];
                if (mymol2 > 1)
                    atomposqq .= myTransData.allAtomPos[:, myAllMolData.cumnumatomslist[mymol2-1]+qq];
                end

                for pp=1:mynumatoms1

                    atompospp .= myTransData.allAtomPos[:, pp];
                    if (mymol1 > 1)
                        atompospp .= myTransData.allAtomPos[:, myAllMolData.cumnumatomslist[mymol1-1]+pp];
                    end
                    atompospp .+= mylatticevec;
                    A[pp, qq] = keffMorse(norm(atompospp .- atomposqq),
                                          k0[pp, qq], a0[pp, qq], Einf[pp, qq]);
                end
            end
        end
    end

end

function assembleDeltaKIdiagblocks!(myDeltaKIsystemdata::DeltaKISystemData{FT},
                                    myAllMolData::MolSystem{OneMol_WithPhonons{FT}}) where FT<:AbstractFloat

    for myDeltaKImolblockdata in myDeltaKIsystemdata.DeltaKImolblockdataarray
    
        mymol1 = myDeltaKImolblockdata.mol1;
        mymol2 = myDeltaKImolblockdata.mol2;

        mynumblocks1 = myDeltaKImolblockdata.numblocks1;
        mynumblocks2 = myDeltaKImolblockdata.numblocks2;
        myblocklist1 = myDeltaKImolblockdata.blocklist1;
        myblocklist2 = myDeltaKImolblockdata.blocklist2;
        
        mynumatoms1 = myAllMolData.numatomslist[mymol1];
        mynumatoms2 = myAllMolData.numatomslist[mymol2];
        
        currDeltaKIdiagelementsblock = view(myDeltaKIsystemdata.DeltaKIdiagelements,
                                            1:myAllMolData.cumnumatomslist[mymol2]);
        if (mymol2 > 1)
            currDeltaKIdiagelementsblock = view(myDeltaKIsystemdata.DeltaKIdiagelements,
                                                myAllMolData.cumnumatomslist[mymol2-1]+1:myAllMolData.cumnumatomslist[mymol2]);
        end

        fill!(currDeltaKIdiagelementsblock, zero(FT));

        for nn2=1:mynumblocks2

            for nn1=1:mynumblocks1

                A = view(myDeltaKImolblockdata.DeltaKImolblock,
                         (nn1-1)*mynumatoms1+1:nn1*mynumatoms1,
                         (nn2-1)*mynumatoms2+1:nn2*mynumatoms2);

                for qq=1:mynumatoms2

                    for pp=1:mynumatoms1

                        currDeltaKIdiagelementsblock[qq] -= A[pp, qq];

                    end

                end

            end

        end

        currDeltaKIdiagelementsblock = view(myDeltaKIsystemdata.DeltaKIdiagelements,
                                            1:myAllMolData.cumnumatomslist[mymol1]);
        if (mymol1 > 1)
            currDeltaKIdiagelementsblock = view(myDeltaKIsystemdata.DeltaKIdiagelements,
                                                myAllMolData.cumnumatomslist[mymol1-1]+1:myAllMolData.cumnumatomslist[mymol1]);
        end

        fill!(currDeltaKIdiagelementsblock, zero(FT));

        for nn2=1:mynumblocks2

            for nn1=1:mynumblocks1

                A = view(myDeltaKImolblockdata.DeltaKImolblock,
                         (nn1-1)*mynumatoms1+1:nn1*mynumatoms1,
                         (nn2-1)*mynumatoms2+1:nn2*mynumatoms2);

                for qq=1:mynumatoms2

                    for pp=1:mynumatoms1

                        currDeltaKIdiagelementsblock[pp] -= A[pp, qq];

                    end

                end

            end

        end

        
    end

end

"""

    assembleDeltaKI!(DeltaKI::AbstractArray{Complex{FT}, 2},
                     myDeltaKIsystemdata::DeltaKISystemData{FT},
                     myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                     myPeriodicData::PeriodicData{FT},
                     blochk::AbstractArray{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the internuclear spring constant matrix of the collection of molecules (which, at
any wavevector `blochk`, will be a Hermitian positive-definite matrix) for the whole system,
including the diagonal contributions, and stamp this into `DeltaKI`.

"""
function assembleDeltaKI!(DeltaKI::AbstractArray{Complex{FT}, 2},
                          myDeltaKIsystemdata::DeltaKISystemData{FT},
                          myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                          myPeriodicData::PeriodicData{FT},
                          blochk::AbstractArray{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    mylatticevec = zeros(FT, 3);
    fill!(DeltaKI, zero(eltype(DeltaKI)));

    # first, do the off-diagonal blocks
    for myDeltaKImolblockdata in myDeltaKIsystemdata.DeltaKImolblockdataarray

        mymol1 = myDeltaKImolblockdata.mol1;
        mymol2 = myDeltaKImolblockdata.mol2;

        if (mymol2 > mymol1) # only assemble the upper-right part, then use Hermiticity for the rest
            startidx1 = 1;
            if (mymol1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[mymol1-1] + 1;
            end
            startidx2 = 1;
            if (mymol2 > 1)
                startidx2 = 3*myAllMolData.cumnumatomslist[mymol2-1] + 1;
            end

            currKIview = view(DeltaKI, startidx1:3*myAllMolData.cumnumatomslist[mymol1],
                              startidx2:3*myAllMolData.cumnumatomslist[mymol2]);
            fill!(currKIview, zero(eltype(currKIview)));

            mynumblocks1 = myDeltaKImolblockdata.numblocks1;
            mynumblocks2 = myDeltaKImolblockdata.numblocks2;
            myblocklist1 = myDeltaKImolblockdata.blocklist1;
            myblocklist2 = myDeltaKImolblockdata.blocklist2;
            
            mynumatoms1 = myAllMolData.numatomslist[mymol1];
            mynumatoms2 = myAllMolData.numatomslist[mymol2];
            
            for nn2=1:mynumblocks2

                for nn1=1:mynumblocks1

                    mylatticevec .= myblocklist1[nn1] .* myPeriodicData.latticevecs[:,1] .+
                    myblocklist2[nn2] .* myPeriodicData.latticevecs[nn2];
                    phaseterm = exp(-1im*dot(blochk, mylatticevec));
                    A = view(myDeltaKImolblockdata.DeltaKImolblock,
                             (nn1-1)*mynumatoms1+1:nn1*mynumatoms1,
                             (nn2-1)*mynumatoms2+1:nn2*mynumatoms2);
                    
                    for qq=1:mynumatoms2

                        for pp=1:mynumatoms1

                            currKIview[3*pp-2:3*pp, 3*qq-2:3*qq] .+=
                            Matrix{Complex{eltype(A)}}(phaseterm * A[pp, qq] * I, 3, 3);

                        end

                    end

                end

            end
            
        end

    end

    # now make it Hermitian, and multiply by 2 to cancel the division by 2
    hermpart!(DeltaKI);
    lmul!(2, DeltaKI);

    # now add the diagonal blocks, which should have no dependence on wavevector
    for qq=1:myAllMolData.cumnumatomslist[end]

        DeltaKI[3*qq-2:3*qq, 3*qq-2:3*qq] .+=
        Matrix{Complex{eltype(myDeltaKIsystemdata.DeltaKIdiagelements)}}(myDeltaKIsystemdata.DeltaKIdiagelements[qq] * I, 3, 3);

    end

    # make it Hermitian again, just to be sure
    hermpart!(DeltaKI);

end
