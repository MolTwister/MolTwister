#
# Copyright (C) 2023 Richard Olsen.
# DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
#
# This file is part of MolTwister.
#
# MolTwister is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MolTwister is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
include(GNUInstallDirs)

# ************************************************************************************
# Function declataions
# ************************************************************************************
function(generate_auto_cuda_files FOLDER_TO_SEARCH GENERATED_LIST_OF_CU_MODULES)
    message(STATUS " ")
    message(STATUS "------------------------------------------------------------------------------------------------------------------------------------------")
    message(STATUS " Automatic generation of CUDA files from CPP files containing the BEGIN_CUDA_COMPATIBLE() macro")
    message(STATUS " Base folder for recursive CPP file search: " ${FOLDER_TO_SEARCH})
    message(STATUS "------------------------------------------------------------------------------------------------------------------------------------------")

    list(APPEND LIST_OF_CU_MODULES ${CMAKE_CURRENT_LIST_DIR}/Utilities/CUDAGeneralizations.h)

    file(GLOB_RECURSE CPP_FILES ${FOLDER_TO_SEARCH}/*.cpp)
    foreach(file ${CPP_FILES})
        message(STATUS ${file})
        file(READ ${file} FILE_CONTENTS)
        string(FIND "${FILE_CONTENTS}" "BEGIN_CUDA_COMPATIBLE()" IS_CUDA_COMPATIBLE)
        if(NOT ${IS_CUDA_COMPATIBLE} EQUAL -1)
            string(REPLACE ".cpp" "_auto_cuda.cu" CUDA_FILE_NAME "${file}")
            string(REPLACE ".cpp" ".h" HEADER_FILE_NAME "${file}")
            message(STATUS "Found BEGIN_CUDA_COMPATIBLE() macro. Hence, the " ${CUDA_FILE_NAME} " file will be created")
            message(STATUS "CUDA module (" ${HEADER_FILE_NAME} ", " ${CUDA_FILE_NAME} ") will be added to the build list for NVCC")
            configure_file(${file} ${CUDA_FILE_NAME} COPYONLY)
            list(APPEND LIST_OF_CU_MODULES ${HEADER_FILE_NAME})
            list(APPEND LIST_OF_CU_MODULES ${CUDA_FILE_NAME})
        endif()
    endforeach(file ${CPP_FILES})

    message(STATUS "------------------------------------------------------------------------------------------------------------------------------------------")
    SET(${GENERATED_LIST_OF_CU_MODULES} ${LIST_OF_CU_MODULES} PARENT_SCOPE)
endfunction(generate_auto_cuda_files)

function(display_list_contents LIST_HEADER LIST_TO_TRAVERS)
    message(STATUS " ")
    message(STATUS "------------------------------------------------------------------------------------------------------------------------------------------")
    message(STATUS " " ${LIST_HEADER})
    message(STATUS "------------------------------------------------------------------------------------------------------------------------------------------")
    set(_LIST_TO_TRAVERS ${LIST_TO_TRAVERS} ${ARGN})
    foreach(item ${_LIST_TO_TRAVERS})
        message(STATUS ${item})
    endforeach()
    message(STATUS "------------------------------------------------------------------------------------------------------------------------------------------")
endfunction()

# ************************************************************************************
# CMake main body
# ************************************************************************************
option(CUDA_DBG_INFO "Include CUDA debug info" OFF)
option(INCLUDE_CUDA_COMMANDS "Include commands that require CUDA" OFF)

if(${INCLUDE_CUDA_COMMANDS})
    project(MolTwister LANGUAGES CXX CUDA)

    # Equivalent to #define INCLUDE_CUDA_COMMANDS=1
    add_definitions(-DINCLUDE_CUDA_COMMANDS=1)
else()
    project(MolTwister LANGUAGES CXX)

    # Equivalent to #define INCLUDE_CUDA_COMMANDS=0
    add_definitions(-DINCLUDE_CUDA_COMMANDS=0)
endif()

set(CMAKE_CXX_STANDARD 17)

find_package(Python3 COMPONENTS Interpreter Development)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Threads REQUIRED)

# Check if Gromacs library is available (used e.g., for XTC trajectory files)
# If so, use it. Note that it can be installed by 'sudo apt install gromacs libgromacs6 libgromacs-dev'
find_library(GROMACS_LIB gromacs NO_CACHE)
if(NOT GROMACS_LIB)
    message(STATUS "Gromacs library not found")
    set(GROMACS_LIB "")
    add_definitions(-DGROMACS_LIB_INSTALLED=0)
else()
    message(STATUS "Gromacs library ${GROMACS_LIB} found")
    add_definitions(-DGROMACS_LIB_INSTALLED=1)
endif()

include_directories(/usr/include/readline)
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${Python3_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
link_libraries(Python3::Python OpenGL::GLU GLUT::GLUT ${GROMACS_LIB})

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
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterTutorialPool.h
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterTutorialPool.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterState.h
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterState.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterVariable.h
    ${CMAKE_CURRENT_LIST_DIR}/MolTwisterVariable.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Utilities/CUDAGeneralizations.h
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
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/XTCFile.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/FFT1D.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/FFT1D.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/FileUtility.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/FileUtility.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/LennardJonesUnits.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/LennardJonesUnits.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/MolDB.h
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/MolDB.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Utilities/Serializer.h

    ${CMAKE_CURRENT_LIST_DIR}/RestAPI/MolTwisterRestAPI.h
    ${CMAKE_CURRENT_LIST_DIR}/RestAPI/MolTwisterRestAPI.cpp

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
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdLl.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdLl.cpp
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
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn.cpp

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
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolTwisterCmdEntry.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/MolTwisterCmdEntry.cpp

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
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdVolumeFromDensity.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdCalculate/CmdVolumeFromDensity.cpp

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
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdFromXTC.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdDcd/CmdFromXTC.cpp

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
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdRadiusOfGyration.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMeasure/CmdRadiusOfGyration.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAtomPos.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAtomPos.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdBondLength.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdBondLength.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdBondType.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdBondType.cpp
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
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdDefaultAtomProps.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdDefaultAtomProps.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAtomName.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdAtomName.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdResname.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdResname.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdUserDefPBC.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdUserDefPBC.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdRotateSel.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMod/CmdRotateSel.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdOptimizeEnergy.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdOptimizeEnergy.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdRun.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdRun.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdFF.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdFF.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdCfg.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdCfg.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/Constants.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/Math.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/Math.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/Particle3D.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/Particle3D.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/NHChain.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/NHChain.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/VelVerlet.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Integrators/VelVerlet.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/MDLoop/MDLoop.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/MDLoop/MDLoop.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/MDLoop/Printf.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/MDLoop/Printf.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/SimulationBox/SimulationBox.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/SimulationBox/SimulationBox.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/ForceFields/MDFFMatrices.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/ForceFields/MDFFMatrices.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/ForceFields/FunctorCalcForce.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/ForceFields/FunctorCalcForce.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/ForceFields/FunctorGenCellList.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/ForceFields/FunctorGenCellList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/ForceFields/FunctorGenNeighList.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/ForceFields/FunctorGenNeighList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Simulator/MDSimulator.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Simulator/MDSimulator.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Config/MolDynConfigStruct.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Config/MolDynConfig.h
    ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/Config/MolDynConfig.cpp
    )

if(${INCLUDE_CUDA_COMMANDS})
    set(MT_CPP_SOURCE_FILES_REQUIRING_CUDA
        ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdCudaTest.h
        ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdCudaTest.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdResetGPU.h
        ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CmdResetGPU.cpp
        )

    generate_auto_cuda_files(${CMAKE_CURRENT_LIST_DIR} AUTO_GENERATED_LIST_OF_CU_MODULES)
    set(MANUALLY_GENERATED_LIST_OF_CU_MODULES
        ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CudaTests/CudaTests.h
        ${CMAKE_CURRENT_LIST_DIR}/Cmd/MolTwisterCmdMolDyn/CudaTests/CudaTests.cu
        ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/CudaDeviceList.h
        ${CMAKE_CURRENT_LIST_DIR}/Cmd/Tools/CudaDeviceList.cu
        )

    set(MT_CUDALIB_SOURCE_FILES ${MANUALLY_GENERATED_LIST_OF_CU_MODULES} ${AUTO_GENERATED_LIST_OF_CU_MODULES})
    display_list_contents("Modules that will be added to the NVCC CUDA build" ${MT_CUDALIB_SOURCE_FILES})
endif()

add_executable(MolTwister ${MT_SOURCE_FILES} ${MT_CPP_SOURCE_FILES_REQUIRING_CUDA})

if(${INCLUDE_CUDA_COMMANDS})
    if(${CUDA_DBG_INFO})
        set(NVIDIA_FLAGS_DBG "${CUDA_NVCC_FLAGS} -g -G -gencode arch=compute_50,code=sm_50")
        set(CUDA_NVCC_FLAGS ${NVIDIA_FLAGS_DBG})
        set(CUDA_NVCC_FLAGS_DEBUG ${NVIDIA_FLAGS_DBG})
        set(CMAKE_CUDA_FLAGS ${NVIDIA_FLAGS_DBG})
        set(CUDA_NVCC_FLAGS_RELEASE ${NVIDIA_FLAGS_DBG})
        message(STATUS "CUDA debug info included...")
    else()
        message(STATUS "CUDA debug info NOT included...")
    endif()

    add_library(MolTwisterCudaLib STATIC ${MT_CUDALIB_SOURCE_FILES})
    target_compile_features(MolTwisterCudaLib PUBLIC cxx_std_14)
    set_target_properties(MolTwisterCudaLib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_property(TARGET MolTwister PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    set(CUDA_MT_LIB MolTwisterCudaLib)
endif()

target_link_libraries(MolTwister ${CUDA_MT_LIB} Threads::Threads Python3::Python GLUT::GLUT OpenGL::GL readline util ${GROMACS_LIB})
target_compile_options(MolTwister PRIVATE -Wno-deprecated -Wno-c++11-extensions -Wno-sign-conversion -Wno-sign-compare -Wno-old-style-cast -Wno-deprecated-declarations -Wno-weak-vtables -Wno-vla-extension -Wno-vla -Wno-padded -Wno-cast-align -Wno-cast-qual -Wno-missing-variable-declarations)

add_custom_target(PythonRestAPI ALL SOURCES "${CMAKE_SOURCE_DIR}/RestAPI/Python/MolTwisterAPI.py")

file(COPY ${CMAKE_SOURCE_DIR}/Tutorials/ DESTINATION ${CMAKE_BINARY_DIR}/moltwister-tutorials/)
file(COPY ${CMAKE_SOURCE_DIR}/RestAPI/Python DESTINATION ${CMAKE_BINARY_DIR}/moltwister-api/)
install(PROGRAMS ${CMAKE_BINARY_DIR}/MolTwister DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} RENAME moltwister)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/Tutorials/ DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/moltwister-tutorials/)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/RestAPI/Python DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBEXECDIR}/moltwister-api/)
