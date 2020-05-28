############################################################################################
module MolModule

using LinearAlgebra;
using DelimitedFiles;
using SharedArrays;

export OneMol, OneMol_WithPhonons, OneMol_NoPhonons, MolSystem, readMolSystem
## export extractKIconst!, extractRequiredMolDataVec, extractOptionalMolDataVec
## export extractRequiredMolDataMat, extractBescale, extractGeneralMolData
## export addKIconsttoKIendatoms!
include("MolData.jl");

end
############################################################################################

############################################################################################
module GenericMatrixVectorUtils

using LinearAlgebra;

include("GenericMatrixVectorUtils.jl");
export addtodiag!, addtodiag, subtractfromdiag!, subtractfromdiag, mulMatDiagleft!
export mulMatDiagleft, mulMatDiagright!, mulMatDiagright, vecNto3N, sympart!, sympart
export hermpart!, hermpart, ahermpart!, ahermpart

end
############################################################################################

############################################################################################
module PeriodicModule

using LinearAlgebra;
using DelimitedFiles;

include("PeriodicData.jl");
export PeriodicData, readPeriodicData, constructlatticefromreciprocal!
export constructreciprocalfromlattice!

end
############################################################################################

############################################################################################
module GFGG

using LinearAlgebra;
using SpecialFunctions;

using ..GenericMatrixVectorUtils;
export addtodiag!, addtodiag, subtractfromdiag!, subtractfromdiag, mulMatDiagleft!
export mulMatDiagleft, mulMatDiagright!, mulMatDiagright, vecNto3N, sympart!, sympart
export hermpart!, hermpart, ahermpart!, ahermpart

include("E1func.jl");
include("GFVACGG.jl");
export c, GFVACGG!, GFVACGG_sca!, GFVACGG0!, GFVACGG0_sca!, GFVACGGcoincident!
export GFVACGG1Ewald!, GFVACGG2Ewald!, GFVACGGEwaldSR, GFVACGGEwaldLR1D, GFVACGGEwaldLR2D
export GFVACGG, GFVACGG_sca, GFVACGG0, GFVACGG0_sca, GFVACGGcoincident, GFVACGG1Ewald
export GFVACGG2Ewald

include("GFPECGG.jl");
export GFPECGG!, GFPECGG0!, GFPECGG1Ewald!, GFPECGG2Ewald!, GFPECGG, GFPECGG0
export GFPECGG1Ewald, GFPECGG2Ewald

end
############################################################################################

############################################################################################
module MolTrans

using LinearAlgebra;
using DelimitedFiles;
using SharedArrays;

using ..MolModule;
export OneMol, OneMol_WithPhonons, OneMol_NoPhonons, MolSystem, readMolSystem

include("MolTransUtils.jl");
export findCOM, translateCOM!, rotatecoords!, translateCOM, rotatecoords, rotate3Nmat!
export rotate3Nmat

include("TransData.jl");
export TransData, readTransData, initializeAtomPos!, transformAtomPos!, untransformAtomPos!

end
############################################################################################

############################################################################################
module MolAlpha

using LinearAlgebra;
using DelimitedFiles;
using SharedArrays;

using ..MolModule;
export OneMol, OneMol_WithPhonons, OneMol_NoPhonons, MolSystem, readMolSystem

using ..GenericMatrixVectorUtils;
export addtodiag!, addtodiag, subtractfromdiag!, subtractfromdiag, mulMatDiagleft!
export mulMatDiagleft, mulMatDiagright!, mulMatDiagright, vecNto3N, sympart!, sympart
export hermpart!, hermpart, ahermpart!, ahermpart

using ..PeriodicModule;
export PeriodicData, readPeriodicData, constructlatticefromreciprocal!
export constructreciprocalfromlattice!

include("MolAlpha.jl");
export constructMolAlpha!, KIk, KIk!

end
############################################################################################

############################################################################################
module DeltaKI

using LinearAlgebra;
using DelimitedFiles;
using SharedArrays;

using ..MolModule;
export OneMol, OneMol_WithPhonons, OneMol_NoPhonons, MolSystem, readMolSystem

using ..GenericMatrixVectorUtils;
export addtodiag!, addtodiag, subtractfromdiag!, subtractfromdiag, mulMatDiagleft!
export mulMatDiagleft, mulMatDiagright!, mulMatDiagright, vecNto3N, sympart!, sympart
export hermpart!, hermpart, ahermpart!, ahermpart

using ..PeriodicModule;
export PeriodicData, readPeriodicData, constructlatticefromreciprocal!
export constructreciprocalfromlattice!

using ..MolTrans;
export findCOM, translateCOM!, rotatecoords!, translateCOM, rotatecoords, rotate3Nmat!
export rotate3Nmat
export TransData, readTransData, initializeAtomPos!, transformAtomPos!, untransformAtomPos!

include("DeltaKIdata.jl");
export DeltaKIMolBlockData, DeltaKISystemData, keffMorse, readDeltaKISystemData
export assembleDeltaKImolblock!, assembleDeltaKIdiagblocks!, assembleDeltaKI!

end
############################################################################################

############################################################################################
module MolEM

using LinearAlgebra;
using DelimitedFiles;
using SharedArrays;
using Distributed;

using ..MolModule;
export OneMol, OneMol_WithPhonons, OneMol_NoPhonons, MolSystem, readMolSystem

include("MolTransUtils.jl");
export findCOM, translateCOM!, rotatecoords!, translateCOM, rotatecoords, rotate3Nmat!
export rotate3Nmat

include("TransData.jl");
export TransData, readTransData, initializeAtomPos!, transformAtomPos!, untransformAtomPos!

using ..GenericMatrixVectorUtils;
export addtodiag!, addtodiag, subtractfromdiag!, subtractfromdiag, mulMatDiagleft!
export mulMatDiagleft, mulMatDiagright!, mulMatDiagright, vecNto3N, sympart!, sympart
export hermpart!, hermpart, ahermpart!, ahermpart

using ..PeriodicModule;
export PeriodicData, readPeriodicData, constructlatticefromreciprocal!
export constructreciprocalfromlattice!

include("MolAlpha.jl");
export constructMolAlpha!, KIk, KIk!

using ..GFGG;
export c, GFVACGG!, GFVACGG_sca!, GFVACGG0!, GFVACGG0_sca!, GFVACGGcoincident!
export GFVACGG1Ewald!, GFVACGG2Ewald!, GFVACGGEwaldSR, GFVACGGEwaldLR1D, GFVACGGEwaldLR2D
export GFVACGG, GFVACGG_sca, GFVACGG0, GFVACGG0_sca, GFVACGGcoincident, GFVACGG1Ewald
export GFVACGG2Ewald
export GFPECGG!, GFPECGG0!, GFPECGG1Ewald!, GFPECGG2Ewald!, GFPECGG, GFPECGG0
export GFPECGG1Ewald, GFPECGG2Ewald

include("MolAssembleGFGG.jl");
export assemble1MolGvacinf!, assemble1MolGvacinfdiag!, subtract1MolalphainvGvac!
export assembleGscaMolBlock!, assembleGvacMolBlock!

include("DeltaKIdata.jl");
export DeltaKIMolBlockData, DeltaKISystemData, keffMorse, readDeltaKISystemData
export assembleDeltaKImolblock!, assembleDeltaKIdiagblocks!, assembleDeltaKI!

end
############################################################################################

############################################################################################
module AARMBEMCCMDARGS

using LinearAlgebra;
using DelimitedFiles;
using Printf;

include("AARMBEMC_ReadCMDARGS.jl");
export readRequiredFilenameARG, readOptionalFilenameARG, readOptionalBoolcondARG
export readOptionalFloatType, readfreqlist!, readklist, readGenvstr

end
############################################################################################

############################################################################################
module AARMBEMCvdW

using LinearAlgebra;
using DelimitedFiles;
using SharedArrays;
using Distributed;
using Printf;

using ..MolEM;
using ..AARMBEMCCMDARGS;

include("vdWenergy.jl");
export vdWenergy

end
############################################################################################

############################################################################################
module AARMBEMCheat

using LinearAlgebra;
using DelimitedFiles;
using SharedArrays;
using Distributed;
using Printf;

using ..MolEM;
using ..AARMBEMCCMDARGS;

include("heattransfer.jl");
export heattransfer

end
############################################################################################
