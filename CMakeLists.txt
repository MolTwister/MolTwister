cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(MolTwister LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)

find_package(PythonLibs REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Threads REQUIRED)

include_directories(/usr/include/readline)
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${PYTHON_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES})

source_group("MolTwister" FILES ${CMAKE_CURRENT_SOURCE_DIR})

set(MT_SOURCE_FILES
    ${CMAKE_CURRENT_LIST_DIR}/main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterPythonExtensions.h
    ${CMAKE_CURRENT_LIST_DIR}/DefaultAtomicProperties.h
    ${CMAKE_CURRENT_LIST_DIR}/DefaultAtomicProperties.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ExpLookup.h
    ${CMAKE_CURRENT_LIST_DIR}/ExpLookup.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MolTwister.h
    ${CMAKE_CURRENT_LIST_DIR}/MolTwister.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterAtom.h
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterAtom.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterCommandPool.h
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterCommandPool.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterState.h
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterState.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterVariable.h
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterVariable.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Utilities/3DVector.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/3DVector.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/3DRect.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/3DRect.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/3DBasis.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/3DBasis.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/ASCIIUtility.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/ASCIIUtility.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/BashColor.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/BashColor.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/ConditionalOnXYZ.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/ConditionalOnXYZ.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/DCDFile.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/DCDFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/FFT1D.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/FFT1D.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/FileUtility.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/FileUtility.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/LennardJonesUnits.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/LennardJonesUnits.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/MolDB.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/MolDB.cpp

    ${CMAKE_CURRENT_LIST_DIR}/3DView/MolTwister3DView.h
    ${CMAKE_CURRENT_LIST_DIR}/3DView/MolTwister3DView.cpp
    ${CMAKE_CURRENT_LIST_DIR}/3DView/MolTwisterGLObject.h
    ${CMAKE_CURRENT_LIST_DIR}/3DView/MolTwisterGLObject.cpp
    ${CMAKE_CURRENT_LIST_DIR}/3DView/MarchingCubes.h
    ${CMAKE_CURRENT_LIST_DIR}/3DView/MarchingCubes.cpp

    ${CMAKE_CURRENT_LIST_DIR}/MDFF/MolTwisterMDFFCoulomb.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/MolTwisterMDFFCoulomb.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/MolTwisterMDFFList.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/MolTwisterMDFFList.cpp

    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBondedList.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBondedList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded_Buck.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded_Buck.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJ.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJ.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJ1208.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJ1208.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJBuck.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJBuck.cpp

    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBond.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBond.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBondList.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBondList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBond_Harm.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBond_Harm.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBond_LJC.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBond_LJC.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBond_Morse.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Bonds/MolTwisterMDFFBond_Morse.cpp

    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Angles/MolTwisterMDFFAngle.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Angles/MolTwisterMDFFAngle.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Angles/MolTwisterMDFFAngleList.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Angles/MolTwisterMDFFAngleList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Angles/MolTwisterMDFFAngle_Class2.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Angles/MolTwisterMDFFAngle_Class2.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Angles/MolTwisterMDFFAngle_Harm.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Angles/MolTwisterMDFFAngle_Harm.cpp

    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Dihedrals/MolTwisterMDFFDih.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Dihedrals/MolTwisterMDFFDih.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Dihedrals/MolTwisterMDFFDihList.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Dihedrals/MolTwisterMDFFDihList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Dihedrals/MolTwisterMDFFDih_Fourier4t.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Dihedrals/MolTwisterMDFFDih_Fourier4t.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Dihedrals/MolTwisterMDFFDih_Harm.h
    ${CMAKE_CURRENT_LIST_DIR}/MDFF/Dihedrals/MolTwisterMDFFDih_Harm.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmd.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmd.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAutoscale.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAutoscale.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCd.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCd.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCopy.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCopy.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDel.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDel.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdGauss9.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdGauss9.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdGenBonds.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdGenBonds.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdGet.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdGet.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdHoomdBlue.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdHoomdBlue.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdLammps.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdLammps.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdList.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdLoad.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdLoad.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdLs.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdLs.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMmol.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMmol.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdPrint.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdPrint.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdPython.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdPython.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdSel.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdSel.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdSet.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdSet.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdVar.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdVar.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdVarlist.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdVarlist.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/ToolsBase.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/ToolsBase.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolTwisterStateTools.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolTwisterStateTools.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolecularTools.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolecularTools.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolecularSystemTools.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolecularSystemTools.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/ParsingTools.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/ParsingTools.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/DCDTools.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/DCDTools.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/ProgressBar.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/ProgressBar.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolTwisterCmdParser.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolTwisterCmdParser.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolTwisterCmdEntry.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolTwisterCmdEntry.h

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdAtom.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdAtom.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdAtoms.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdAtoms.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdMDNonBonded.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdMDNonBonded.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdMDBond.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdMDBond.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdMDAngle.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdMDAngle.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdMDDihedral.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdAdd/CmdMDDihedral.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdQBal.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdQBal.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDensityProfile.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDensityProfile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdVACF.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdVACF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdFFT.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdFFT.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdVDOS.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdVDOS.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdMSD.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdMSD.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdPairCorrelation.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdPairCorrelation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdCOM.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdCOM.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdPotenEnergyMap.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdPotenEnergyMap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdForceBetween.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdForceBetween.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdEnergyBetween.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdEnergyBetween.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDensityMap.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDensityMap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdLoading.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdLoading.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdEnergyOfTranslation.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdEnergyOfTranslation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdHBondCount.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdHBondCount.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDihedralDistr.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDihedralDistr.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDihedralDistrCOM.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDihedralDistrCOM.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDistProbabilityCOM.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDistProbabilityCOM.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMomProfile.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMomProfile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMomPerturbCharge.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMomPerturbCharge.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMomPerturbZPosExchange.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMomPerturbZPosExchange.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMomPerturbZPos.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMomPerturbZPos.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMom.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdDipMom.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdReadRecord.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdReadRecord.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdNumRecords.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdNumRecords.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdReadCoordinate.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdReadCoordinate.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdNumCoordinates.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdNumCoordinates.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdHeader.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdHeader.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdUnwrap.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdUnwrap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdAtomicUnwrap.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdAtomicUnwrap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdWrap.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdWrap.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdBondLength.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdBondLength.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdBondLengthDyn.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdBondLengthDyn.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdBondSep.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdBondSep.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdAngle.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdAngle.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdDihedral.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdDihedral.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCoulombEnergy.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCoulombEnergy.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCoulombPotential.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCoulombPotential.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdAtomPos.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdAtomPos.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCenter.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCenter.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdPBC.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdPBC.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCharge.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCharge.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCount.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdCount.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdBondCount.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdBondCount.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdOverlap.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdOverlap.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAtomPos.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAtomPos.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdBondLength.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdBondLength.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAngle.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAngle.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdDihedral.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdDihedral.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdCharge.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdCharge.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdMass.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdMass.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdMobillity.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdMobillity.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdSigma.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdSigma.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAtomName.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAtomName.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdResname.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdResname.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdUserDefPBC.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdUserDefPBC.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdRotateSel.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdRotateSel.cpp
    )

add_executable(MolTwister ${MT_SOURCE_FILES})
target_link_libraries(MolTwister Threads::Threads ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} readline)
target_compile_options(MolTwister PRIVATE -Wno-deprecated -Wno-c++11-extensions -Wno-sign-conversion -Wno-sign-compare -Wno-old-style-cast -Wno-deprecated-declarations -Wno-weak-vtables -Wno-vla-extension -Wno-vla -Wno-padded -Wno-cast-align -Wno-cast-qual -Wno-missing-variable-declarations)
