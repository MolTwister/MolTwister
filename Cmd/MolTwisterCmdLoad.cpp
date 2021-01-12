//
// Copyright (C) 2021 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include "Utilities/FileUtility.h"
#include "MolTwisterCommandPool.h"
#include "MolTwisterCmdLoad.h"
#include "Tools/MolTwisterStateTools.h"
#include "Tools/ProgressBar.h"
#include "../Utilities/DCDFile.h"

void CCmdLoad::onAddKeywords()
{
    addKeyword("load");
    addKeyword("xyz");
    addKeyword("dcd");
    addKeyword("pdb");
    addKeyword("mtt");
    addKeyword("script");
    addKeyword("masscharge");
    addKeyword("qepos");
    addKeyword("ignore");
    addKeyword("all");
    addKeyword("frame");
    addKeyword("noquery");
}

std::string CCmdLoad::getHelpString() const
{ 
    std::string  szText;
    
    szText+= "\tUsage: load <data type> <file path>\r\n";
    szText+= "\r\n";
    szText+= "\tLoad data into MolTwister. The allowed data types are:\r\n";
    szText+= "\r\n";
    szText+= "\txyz :        XYZ-coordinate files (including XYZ files containing\r\n";
    szText+= "\t             several frames). Bonds can be ignored using 'ignore' (\r\n";
    szText+= "\t             see genbonds). Use 'frame <frame>' to select frame to use\r\n";
    szText+= "\t             as basis for generating the bonds. Default is frame zero.\r\n";
    szText+= "\t             The order of the additional arguments is 'ignore', 'frame'.\r\n";
    szText+= "\tdcd :        DCD-coordinate files. Bonds can be ignored using 'ignore' (\r\n";
    szText+= "\t             see genbonds). Use 'frame <frame>' to select frame to use\r\n";
    szText+= "\t             as basis for generating the bonds. Default is frame zero.\r\n";
    szText+= "\t             Use 'stride <stride>' to only load every <stride> frame.\r\n";
    szText+= "\t             The order of the additional arguments is 'ignore', 'frame', 'stride'.\r\n";
    szText+= "\tpdb :        PDB-coordinate files. Bonds can be ignored using 'ignore'.\r\n";
    szText+= "\t             To avoid query about delete ('del'), append ('app') or cancel ('can'),\r\n";
    szText+= "\t             use 'noquery <response>' with the appropriate response.\r\n";
    szText+= "\tmtt :        MolTwister trajectory file.\r\n";
    szText+= "\tscript :     MolTwister script containing a sequence of MolTwister commands,\r\n";
    szText+= "\t             each formated as: exec(\"<MolTwister command>\");. For example,\r\n";
    szText+= "\t             exec(\"add atom C at 0 0 0\"); would add a C atom at (0, 0, 0)\r\n";
    szText+= "\tmasscharge : Load charges and masses into the current molecular structure. In\r\n";
    szText+= "\t             the file each line defines a charge/mass in two possible ways:\r\n";
    szText+= "\t                   * ID <ID, e.g. H, C, C2, or similar> <partial charge> <mass>\r\n";
    szText+= "\t                   * AtInd <Index of atom> <partial charge> <mass>\r\n";
    szText+= "\tqepos :      XYZ-position files (*.pos) from QuantumEspresso. Note that in\r\n";
    szText+= "\t             this case <file path> = <*.pos file> <input file>\r\n\r\n";
    szText+= "\t             If no data type is selected MolTwister will try to recognize\r\n";
    szText+= "\t             the file format from the file extension.\r\n";
    szText+= "\tpython :     Execute a python script. How to access MolTwister from a Python\r\n";
    szText+= "\t             is documented under the 'mtpython' command. However, when loading\r\n";
    szText+= "\t             a script, only the contents of 'mtpython {<contents>}' is to be\r\n";
    szText+= "\t             included in the Python script.";

    return szText; 
}

void CCmdLoad::execute(std::string commandLine)
{
    std::vector<std::string> bondAtomsToIgnore;
    bool genBonds = true;
    bool ignoreAllBonds = false;
    int baseFrameIndex = 0;
    int arg = 1;
    
    if(!state_) return;
    
    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    if((text != "xyz") &&
       (text != "dcd") &&
       (text != "pdb") &&
       (text != "mtt") &&
       (text != "script") &&
       (text != "python"))
    {
        // Try and see if szText1 is a recognized filename
        // and not a specification of file type
        std::string ext = CFileUtility::getExtension(text);
        if(ext == "xyz") { text = "xyz"; arg--; }
        if(ext == "dcd") { text = "dcd"; arg--; }
        if(ext == "pdb") { text = "pdb"; arg--; }
        if(ext == "mtt") { text = "mtt"; arg--; }
        if(ext == "script") { text = "script"; arg--; }
        if(ext == "py") { text = "python"; arg--; }
    }
    
    if(text == "xyz")
    {
        parseXYZCommand(commandLine, arg, bondAtomsToIgnore, genBonds, ignoreAllBonds, baseFrameIndex);
    }
    else if(text == "dcd")
    {
        parseDCDCommand(commandLine, arg, bondAtomsToIgnore, genBonds, ignoreAllBonds, baseFrameIndex);
    }
    else if(text == "pdb")
    {
        parsePDBCommand(commandLine, arg, bondAtomsToIgnore, genBonds, ignoreAllBonds, baseFrameIndex);
    }
    else if(text == "mtt")
    {
        parseMTTCommand(commandLine, arg, bondAtomsToIgnore, genBonds, ignoreAllBonds, baseFrameIndex);
    }
    else if(text == "script")
    {
        parseScriptCommand(commandLine, arg, bondAtomsToIgnore, genBonds, ignoreAllBonds, baseFrameIndex);
    }
    else if(text == "python")
    {
        parsePythonCommand(commandLine, arg, bondAtomsToIgnore, genBonds, ignoreAllBonds, baseFrameIndex);
    }
    else if(text == "masscharge")
    {
        parseMasschargeCommand(commandLine, arg, bondAtomsToIgnore, genBonds, ignoreAllBonds, baseFrameIndex);
    }
    else if(text == "qepos")
    {
        parseQeposCommand(commandLine, arg, bondAtomsToIgnore, genBonds, ignoreAllBonds, baseFrameIndex);
    }
    else
    {
        printf("Syntax Error: First argument should be the file type!");
    }
    
    if(genBonds && (state_->atoms_.size() > 0))
    {
        if(ignoreAllBonds)
            state_->searchForAtomTypes(bondAtomsToIgnore);
        if((baseFrameIndex < 0) || (baseFrameIndex >= (int)state_->atoms_[0]->r_.size()))
        {
            for(int i=0; i<(int)state_->atoms_[0]->r_.size(); i++)
                CMolTwisterStateTools(state_, stdOut_).generateBonds(0.8, false, true, i, nullptr, &bondAtomsToIgnore);
        }
        else
        {
            CMolTwisterStateTools(state_, stdOut_).generateBonds(0.8, false, true, baseFrameIndex, nullptr, &bondAtomsToIgnore);
        }
    }
    
    if(state_->view3D_) state_->view3D_->requestUpdate(true);    
}

void CCmdLoad::parseXYZCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool&, int& baseFrameIndex)
{
    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    if(text.length() == 0)
    {
        printf("Syntax Error: Second argument should be the filename!");
    }
    else
    {
        std::string ignoreString = CASCIIUtility::getWord(commandLine, arg++);
        CASCIIUtility::removeWhiteSpace(ignoreString);
        if(ignoreString == "ignore")
        {
            ignoreString = CASCIIUtility::getWord(commandLine, arg++);
            CASCIIUtility::removeWhiteSpace(ignoreString);
            bondAtomsToIgnore = CASCIIUtility::getWords(ignoreString, ",");
        }
        else arg--;

        std::string baseFrame;
        CASCIIUtility::removeWhiteSpace(baseFrame);
        if(baseFrame == "frame")
        {
            baseFrame = CASCIIUtility::getWord(commandLine, arg++);
            CASCIIUtility::removeWhiteSpace(baseFrame);
            baseFrameIndex = atoi(baseFrame.data());
        }
        else
        {
            baseFrameIndex = 0;
            arg--;
        }

        if(text[0] == '/')
        {
            if(!readXYZFile(text.data(), genBonds))
                printf("Error loading XYZ file!");
        }
        else
        {
            std::string filePath = CFileUtility::getCWD() + text;

            if(!readXYZFile(filePath.data(), genBonds))
                printf("Error loading XYZ file!");
        }
    }
}

void CCmdLoad::parseDCDCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool&, int& baseFrameIndex)
{
    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    if(text.length() == 0)
    {
        printf("Syntax Error: Second argument should be the filename!");
    }
    else
    {
        std::string ignoreString = CASCIIUtility::getWord(commandLine, arg++);
        CASCIIUtility::removeWhiteSpace(ignoreString);
        if(ignoreString == "ignore")
        {
            ignoreString = CASCIIUtility::getWord(commandLine, arg++);
            CASCIIUtility::removeWhiteSpace(ignoreString);
            bondAtomsToIgnore = CASCIIUtility::getWords(ignoreString, ",");
        }
        else arg--;

        std::string baseFrame;
        CASCIIUtility::removeWhiteSpace(baseFrame);
        if(baseFrame == "frame")
        {
            baseFrame = CASCIIUtility::getWord(commandLine, arg++);
            CASCIIUtility::removeWhiteSpace(baseFrame);
            baseFrameIndex = atoi(baseFrame.data());
        }
        else
        {
            baseFrameIndex = 0;
            arg--;
        }

        int stride = 1;
        std::string strideString = CASCIIUtility::getWord(commandLine, arg++);
        CASCIIUtility::removeWhiteSpace(strideString);
        if(strideString == "stride")
        {
            stride = atoi(CASCIIUtility::getWord(commandLine, arg++).data());
        }
        else arg++;

        if(text[0] == '/')
        {
            if(!readDCDFile(text.data(), genBonds, stride))
                printf("Error loading DCD file!");
        }
        else
        {
            std::string filePath = CFileUtility::getCWD() + text;

            if(!readDCDFile(filePath.data(), genBonds, stride))
                printf("Error loadingDCD file!");
        }
    }
}

void CCmdLoad::parsePDBCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int&)
{
    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    if(text.length() == 0)
    {
        printf("Syntax Error: Second argument should be the filename!");
    }
    else
    {
        std::string typeString;

        typeString = CASCIIUtility::getWord(commandLine, arg++);
        CASCIIUtility::removeWhiteSpace(typeString);
        if(typeString == "ignore")
        {
            typeString = CASCIIUtility::getWord(commandLine, arg++);
            if(typeString == "all")
            {
                ignoreAllBonds = true;
            }
            else
            {
                CASCIIUtility::removeWhiteSpace(typeString);
                bondAtomsToIgnore = CASCIIUtility::getWords(typeString, ",");
            }
        }
        else arg--;

        typeString = CASCIIUtility::getWord(commandLine, arg++);
        CASCIIUtility::removeWhiteSpace(typeString);
        auto noQuery = std::pair<bool, std::string>(false, "");
        if(typeString == "noquery")
        {
            noQuery.first = true;
            noQuery.second = CASCIIUtility::getWord(commandLine, arg++);
        }
        else arg--;

        if(text[0] == '/')
        {
            if(!readPDBFile(text.data(), genBonds, noQuery))
                printf("Error loading PDB file!");
        }
        else
        {
            std::string filePath = CFileUtility::getCWD() + text;

            if(!readPDBFile(filePath.data(), genBonds, noQuery))
                printf("Error loading PDB file!");
        }
    }
}

void CCmdLoad::parseMTTCommand(std::string commandLine, int& arg, std::vector<std::string>&, bool& genBonds, bool&, int&)
{
    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    if(text.length() == 0)
    {
        printf("Syntax Error: Second argument should be the filename!");
    }
    else
    {
        if(text[0] == '/')
        {
            if(!readMTTFile(text.data()))
                printf("Error loading MTT file!");
        }
        else
        {
            std::string filePath = CFileUtility::getCWD() + text;

            if(!readMTTFile(filePath.data()))
                printf("Error loading MTT file!");
        }
    }

    genBonds = false;
}

void CCmdLoad::parseScriptCommand(std::string commandLine, int& arg, std::vector<std::string>&, bool& genBonds, bool&, int&)
{
    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    if(text.length() == 0)
    {
        printf("Syntax Error: Second argument should be the filename!");
    }
    else
    {
        if(text[0] == '/')
        {
            if(!readScriptFile(text.data()))
                printf("Error loading Script file!");
        }
        else
        {
            std::string filePath = CFileUtility::getCWD() + text;

            if(!readScriptFile(filePath.data()))
                printf("Error loading Script file!");
        }
    }

    genBonds = false;
}

void CCmdLoad::parsePythonCommand(std::string commandLine, int& arg, std::vector<std::string>&, bool& genBonds, bool&, int&)
{
    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    if(text.length() == 0)
    {
        printf("Syntax Error: Second argument should be the filename!");
    }
    else
    {
        if(text[0] == '/')
        {
            if(!readPythonFile(text.data()))
                printf("Error loading Python file!");
        }
        else
        {
            std::string filePath = CFileUtility::getCWD() + text;

            if(!readPythonFile(filePath.data()))
                printf("Error loading Python file!");
        }
    }

    genBonds = false;
}

void CCmdLoad::parseMasschargeCommand(std::string commandLine, int& arg, std::vector<std::string>&, bool& genBonds, bool&, int&)
{
    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    if(text.length() == 0)
    {
        printf("Syntax Error: Second argument should be the filename!");
    }
    else
    {
        if(text[0] == '/')
        {
            if(!readMassChargeFile(text.data()))
                printf("Error loading Mass-Charge file!");
        }
        else
        {
            std::string filePath = CFileUtility::getCWD() + text;

            if(!readMassChargeFile(filePath.data()))
                printf("Error loading Mass-Charge file!");

        }
    }

    genBonds = false;
}

void CCmdLoad::parseQeposCommand(std::string commandLine, int& arg, std::vector<std::string>&, bool&, bool&, int&)
{
    std::string text1 = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text1);
    std::string text2 = CASCIIUtility::getWord(commandLine, arg++);
    CASCIIUtility::removeWhiteSpace(text2);
    if((text1.length() == 0) || (text2.length() == 0))
    {
        printf("Syntax Error: Second argument should be the *.pos and third should be the input filename!");
    }
    else
    {
        if(text1[0] != '/')
        {
            std::string filePath = CFileUtility::getCWD() + text1;

            text1 = filePath;
        }

        if(text2[0] != '/')
        {
            std::string filePath = CFileUtility::getCWD() + text1;

            text2 = filePath;
        }

        if(!readQEPosFile(text1.data(), text2.data()))
            printf("Error loading QuantumEspresso *.pos or input file!");
    }
}

bool CCmdLoad::readXYZFile(std::string xyzFileName, bool& genBonds)
{
    CProgressBar progBar;
    double X, Y, Z;
    std::string ID, line, text;
    FILE* fileHandle = fopen(xyzFileName.data(), "r");
    bool lastLineInFile, firstFrameAdded=false;
    int numCoordinates;

    
    if(!state_) return false;
    
    if(!fileHandle)
    {
        printf("\r\nCould not open file: %s!\r\n", xyzFileName.data());
        return false;
    }

    
    // Get confiramation from user and purge everything if 'Yes'
    genBonds = false;
    printf("Delete all atoms (del), Append (app) or Cancel (can)?: ");
    char stringAnswer[20];
    fgets(stringAnswer, 19, stdin);
    std::string action = CASCIIUtility::removeCRLF(stringAnswer);
    if(action == "can")        return true;
    else if(action == "app")   { /*no action required*/ }
    else if(action == "del")   state_->purgeAtomsList();
    else                       return true;


    // Count lines in the file
    int numLinesInFile = 0;
    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        numLinesInFile++;

    } while(!lastLineInFile);
    fseek(fileHandle, 0, SEEK_SET);

    
    // Read coordinates from file
    progBar.beginProgress("Loading XYZ file contents");
    int lineIndex = 0;
    do
    {
        // Read num coordinates in XYZ file
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        lineIndex++;

        CASCIIUtility::removeWhiteSpace(line);
        numCoordinates = atoi(line.data());

        // Read comment line and positions of XYZ file
        if(!lastLineInFile && (numCoordinates > 0))
        {
            int indexFrame=0;

            if(firstFrameAdded) indexFrame = state_->addFrame();
            
            // Read comment line of XYZ file
            line = CFileUtility::readLine(fileHandle, lastLineInFile);
            lineIndex++;
            
            // Read positions
            for(int i=0; i<numCoordinates; i++)
            {
                line = CFileUtility::readLine(fileHandle, lastLineInFile);
                lineIndex++;
                
                ID = CASCIIUtility::getWord(line, 0);
                
                text = CASCIIUtility::getWord(line, 1);
                X = atof(text.data());
                
                text = CASCIIUtility::getWord(line, 2);
                Y = atof(text.data());
                
                text = CASCIIUtility::getWord(line, 3);
                Z = atof(text.data());
                
                if(!firstFrameAdded)
                    state_->addAtom(X, Y, Z, ID.data());
                else
                {
                    if(indexFrame >= 0) state_->setAtomCoordinates(indexFrame, i, X, Y, Z);
                }
                if(i < (int)state_->atoms_.size())
                    state_->atoms_[i]->sigma_ = state_->defaultAtProp_.getWDWRadius(ID.data());
            }

            firstFrameAdded = true;
        }

        progBar.updateProgress(lineIndex, numLinesInFile);
   
    } while(!lastLineInFile);
    progBar.endProgress();

    genBonds = true;
    fclose(fileHandle);
    return true;
}

bool CCmdLoad::readDCDFile(std::string dcdFileName, bool& genBonds, int stride)
{
    CProgressBar progBar;
    double X, Y, Z;
    std::string ID;

    genBonds = false;
    CDCDFile dcdFile;

    if(!dcdFile.open(dcdFileName.data(), stride))
    {
        printf("\r\nCould not open file: %s!\r\n", dcdFileName.data());
        return false;
    }

    int numAtoms = dcdFile.getNumCoordinatesInRecord();
    if((int)state_->atoms_.size() !=  numAtoms)
    {
        printf("\r\nThe file, %s! It was found to have %i atoms per record (based on first record in file), while the loaded system has %i atoms. The number of atoms much match!\r\n", dcdFileName.data(), numAtoms, (int)state_->atoms_.size());
        return false;
    }

    progBar.beginProgress("Loading DCD file contents");
    int numRecords = dcdFile.getNumRecords();
    bool firstFrameAdded = false;
    for(int i=0; i<numRecords; i++)
    {
        dcdFile.gotoRecord(i);

        int indexFrame=0;
        if(firstFrameAdded) indexFrame = state_->addFrame();

        numAtoms = dcdFile.getNumCoordinatesInRecord();
        const float* coordDataX = dcdFile.getCoordinateDataX();
        const float* coordDataY = dcdFile.getCoordinateDataY();
        const float* coordDataZ = dcdFile.getCoordinateDataZ();
        for(int j=0; j<numAtoms; j++)
        {
            X = (double)coordDataX[j];
            Y = (double)coordDataY[j];
            Z = (double)coordDataZ[j];

            ID = state_->atoms_[j]->getID();

            if(indexFrame >= 0) state_->setAtomCoordinates(indexFrame, j, X, Y, Z);
            if(j < (int)state_->atoms_.size())
                state_->atoms_[j]->sigma_ = state_->defaultAtProp_.getWDWRadius(ID.data());
        }

        firstFrameAdded = true;
        progBar.updateProgress(i+1, numRecords);
    }
    progBar.endProgress();

    genBonds = true;
    return true;
}

bool CCmdLoad::readPDBFile(std::string pdbFileName, bool& genBonds, const std::pair<bool, std::string>& noQuery)
{
    std::string stringPart;
    double X, Y, Z;
    std::string ID, resname, line, text;
    FILE* fileHandle = fopen(pdbFileName.data(), "r");
    bool lastLineInFile, foundEnd, isHetAtm, isAtom;
    int  coordIndex = 0;
 
    
    if(!state_) return false;
    
    if(!fileHandle)
    {
        printf("\r\nCould not open file: %s!\r\n", pdbFileName.data());
        return false;
    }
    
    
    // Get confiramation from user and purge everything if 'del'
    genBonds = false;
    std::string action;
    if(!noQuery.first)
    {
        printf("Delete all atoms (del), Append (app) or Cancel (can)?: ");
        char stringAnswer[20];
        fgets(stringAnswer, 19, stdin);
        action = CASCIIUtility::removeCRLF(stringAnswer);
    }
    else
    {
        printf("Delete all atoms (del), Append (app) or Cancel (can)?: %s\r\n", noQuery.second.data());
        action = noQuery.second;
    }
    if(action == "can")        return true;
    else if(action == "app")   { /*no action required*/ }
    else if(action == "del")   state_->purgeAtomsList();
    else                       return true;


    // Read coordinates from file
    do
    {
        // Read positions
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        foundEnd = (CASCIIUtility::findString("END", line) != -1) ? true : false;
        if(!foundEnd)
        {
            stringPart = CASCIIUtility::extractString(0, 5, line);
            CASCIIUtility::removeWhiteSpace(stringPart);
            isHetAtm = (stringPart == "HETATM") ? true : false;
            isAtom = (stringPart == "ATOM") ? true : false;
            
            if(isHetAtm || isAtom)
            {
                stringPart = CASCIIUtility::extractString(12, 15, line);
                CASCIIUtility::removeWhiteSpace(stringPart);
                ID = stringPart;

                stringPart = CASCIIUtility::extractString(17, 19, line);
                CASCIIUtility::removeWhiteSpace(stringPart);
                resname = stringPart;
                
                stringPart = CASCIIUtility::extractString(30, 37, line);
                X = atof(stringPart.data());
                
                stringPart = CASCIIUtility::extractString(38, 45, line);
                Y = atof(stringPart.data());
                
                stringPart = CASCIIUtility::extractString(46, 53, line);
                Z = atof(stringPart.data());
                
                coordIndex = state_->addAtom(X, Y, Z, ID.data());
                
                if(coordIndex < (int)state_->atoms_.size())
                {
                    state_->atoms_[coordIndex]->sigma_ = state_->defaultAtProp_.getWDWRadius(ID.data());
                    state_->atoms_[coordIndex]->resname_ = resname;
                }
            }
        }
        
    } while(!lastLineInFile && !foundEnd);

    genBonds = true;
    fclose(fileHandle);
    return true;
}

bool CCmdLoad::readMTTFile(std::string mttFileName)
{
    std::string line, text;
    FILE* fileHandle = fopen(mttFileName.data(), "r");
    bool lastLineInFile;
    
    
    if(!state_) return false;
    
    if(!fileHandle)
    {
        printf("\r\nCould not open file: %s!\r\n", mttFileName.data());
        return false;
    }
    
    
    // Get confiramation from user and purge everything if 'del'
    printf("Delete all atoms (del), Append (app) or Cancel (can)?: ");
    char stringAnswer[20];
    fgets(stringAnswer, 19, stdin);
    std::string action = CASCIIUtility::removeCRLF(stringAnswer);
    if(action == "can")        return true;
    else if(action == "app")   { /*no action required*/ }
    else if(action == "del")   state_->purgeAtomsList();
    else                       return true;

    
    // Read file header
    line = CFileUtility::readLine(fileHandle, lastLineInFile);
    text = CASCIIUtility::getWord(line, 0);
    if(text != "mtt")
    {
        printf("Error! %s mot a MolTwister trajectory (MTT) file!", mttFileName.data());
        return true;
    }

    line = CFileUtility::readLine(fileHandle, lastLineInFile);
    text = CASCIIUtility::getWord(line, 0);
    unsigned long ulVer = atoi(text.data());
    if(ulVer > 1)
    {
        printf("Error! this version of MolTwister does not support version %lu MTT files (version must be <= 1)!", ulVer);
        return true;
    }

    line = CFileUtility::readLine(fileHandle, lastLineInFile);
    text = CASCIIUtility::getWord(line, 0);
    bool isBinary = (text == "binary") ? true : false;
    
    bool molIndicesAvailable;
    bool resnamesAvailable;
    bool bondsAvailable;
    if(!isBinary)
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        text = CASCIIUtility::getWord(line, 0);
        molIndicesAvailable = (atoi(text.data()) == 0) ? false : true;

        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        text = CASCIIUtility::getWord(line, 0);
        resnamesAvailable = (atoi(text.data()) == 0) ? false : true;
        
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        text = CASCIIUtility::getWord(line, 0);
        bondsAvailable = (atoi(text.data()) == 0) ? false : true;
    }
    else
    {
        unsigned char byteRead;
        fread(&byteRead, sizeof(byteRead), 1, fileHandle);
        molIndicesAvailable = (byteRead == 0) ? false : true;
        fread(&byteRead, sizeof(byteRead), 1, fileHandle);
        resnamesAvailable = (byteRead == 0) ? false : true;
        fread(&byteRead, sizeof(byteRead), 1, fileHandle);
        bondsAvailable = (byteRead == 0) ? false : true;
    }
    

    // Read start delimiter of coordinate frame
    int numPrevLoadMolecules = CMolTwisterStateTools(state_, stdOut_).getNumMolecules();
    int numPrevLoadAtoms = (int)state_->atoms_.size();
    bool foundOneMoreRecord = false;
    if(!isBinary)
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        text = CASCIIUtility::getWord(line, 0);
        if(text == "{") foundOneMoreRecord = true;
    }
    else
    {
        long longRead;
        fread(&longRead, sizeof(longRead), 1, fileHandle);
        if(longRead == -1) foundOneMoreRecord = true;
    }
    
    if(foundOneMoreRecord)
    {
        // Read coordinate frame. Stop at end delimiter, or end of file
        int coordIndex = 0;
        bool end = false;
        do
        {
            long molIndex = 0;
            double X = 0.0, Y = 0.0, Z = 0.0;
            std::string ID, sresname = "-";
            std::vector<long> bondsTo;
            if(!isBinary)
            {
                int arg = 0;
                line = CFileUtility::readLine(fileHandle, lastLineInFile);
                text = CASCIIUtility::getWord(line, arg++);
                if(text == "}") end = true;
                if(!end)
                {
                    text = CASCIIUtility::getWord(line, arg++);
                    if(text == "---") ID = "";
                    else              ID = text;

                    text = CASCIIUtility::getWord(line, arg++);
                    X = atof(text.data());
                    text = CASCIIUtility::getWord(line, arg++);
                    Y = atof(text.data());
                    text = CASCIIUtility::getWord(line, arg++);
                    Z = atof(text.data());
                    
                    if(molIndicesAvailable)
                    {
                        text = CASCIIUtility::getWord(line, arg++);
                        molIndex = atoi(text.data());
                    }
                    
                    if(resnamesAvailable)
                    {
                        text = CASCIIUtility::getWord(line, arg++);
                        if(text == "---")  sresname = "";
                        else               sresname = text;
                    }
                    
                    if(bondsAvailable)
                    {
                        text = CASCIIUtility::getWord(line, arg++);
                        unsigned char byteNumBonds = (unsigned char)atoi(text.data());
                        
                        for(unsigned char l=0; l<byteNumBonds; l++)
                        {
                            text = CASCIIUtility::getWord(line, arg++);
                            bondsTo.emplace_back(atoi(text.data()));
                        }
                    }
                }
                
                if(lastLineInFile) end = true;
            }
            else
            {
                long longRead;
                unsigned char byteRead;
                if(fread(&longRead, sizeof(longRead), 1, fileHandle) != 1) end = true;
                if(longRead == -2) end = true;
                if(!end)
                {
                    fread(&byteRead, sizeof(byteRead), 1, fileHandle);
                    for(int i=0; i<byteRead; i++)
                    {
                        char character;
                        fread(&character, sizeof(character), 1, fileHandle);
                        ID+= character;
                    }
                    
                    fread(&X, sizeof(X), 1, fileHandle);
                    fread(&Y, sizeof(Y), 1, fileHandle);
                    fread(&Z, sizeof(Z), 1, fileHandle);
                    
                    if(molIndicesAvailable)
                        fread(&molIndex, sizeof(molIndex), 1, fileHandle);
                    
                    if(resnamesAvailable)
                    {
                        fread(&byteRead, sizeof(byteRead), 1, fileHandle);
                        for(int i=0; i<byteRead; i++)
                        {
                            char character;
                            fread(&character, sizeof(character), 1, fileHandle);
                            sresname+= character;
                        }
                    }
                    
                    if(bondsAvailable)
                    {
                        unsigned char numBonds;
                        fread(&numBonds, sizeof(numBonds), 1, fileHandle);
                        for(unsigned char l=0; l<numBonds; l++)
                        {
                            fread(&longRead, sizeof(longRead), 1, fileHandle);
                            bondsTo.emplace_back(longRead);
                        }
                    }
                }
            }
            
            if(!end)
            {
                coordIndex = state_->addAtom(X, Y, Z, ID.data());
                
                if(coordIndex < (int)state_->atoms_.size())
                {
                    state_->atoms_[coordIndex]->sigma_ = state_->defaultAtProp_.getWDWRadius(ID.data());
                    state_->atoms_[coordIndex]->resname_ = sresname;
                    
                    if(molIndicesAvailable)
                    {
                        state_->atoms_[coordIndex]->setMolIndex((int)molIndex + numPrevLoadMolecules);
                    }
                    
                    if(bondsAvailable)
                    {
                        for(int l=0; l<(int)bondsTo.size(); l++)
                        {
                            int atomToConnectTo = numPrevLoadAtoms + (int)bondsTo[l];
                            if(atomToConnectTo < (int)state_->atoms_.size())
                            {
                                CAtom* pAtom = state_->atoms_[atomToConnectTo].get();
                                state_->atoms_[coordIndex]->attachBondTo(pAtom);
                                pAtom->attachBondTo(state_->atoms_[coordIndex].get());
                            }
                        }
                    }
                }
            }

        } while(!end);
    }
    
    fclose(fileHandle);
    return true;
}

bool CCmdLoad::readScriptFile(std::string scriptFileName)
{
    bool lastLineInFile = false;
    FILE* fileHandle = fopen(scriptFileName.data(), "r");
    std::string line, command, funcName, argument;
    std::vector<std::shared_ptr<CCmd>> cmdList;

    if(!state_) return false;

    if(!fileHandle)
    {
        printf("\r\nError, could not open file: %s!\r\n", scriptFileName.data());
        return false;
    }
    
    CMolTwisterCommandPool::generateCmdList(state_, cmdList);

    fprintf(stdOut_, "\r\n");
    
    do
    {
        line = CFileUtility::readToNextDelimiterIgnoreCppComments(fileHandle, ';', lastLineInFile);
        CASCIIUtility::parseCLikeFuncCall(line, funcName, argument);
    
        if(funcName == "exec")
        {
            CASCIIUtility::removeWhiteSpace(argument, "\"");
            command = CASCIIUtility::getWord(argument, 0);

            int pipeSymbIndex = CASCIIUtility::findString("> ", argument);
            
            for(int i=0; i<(int)cmdList.size(); i++)
            {
                if(cmdList[i] && cmdList[i]->checkCmd(command.data()))
                {
                    std::string fileName;

                    // Check if '> <file>' (i.e. pipe) was requested. If so
                    // redirect all fprintf() output to <file> and keep all printf()
                    // statements directed to stdout
                    FILE* stdOutFile = stdout;
                    if(pipeSymbIndex != -1)
                    {
                        fileName = argument.substr(pipeSymbIndex+1, std::string::npos);
                        CASCIIUtility::removeWhiteSpace(fileName);
                        stdOutFile = fopen(fileName.data(), "w+");
                        
                        if(stdOutFile)  cmdList[i]->redirectOutput(stdOutFile);
                        else            printf("Error: could not create file %s!", fileName.data());
                    }
                    
                    // Execute command
                    fprintf(stdOut_, "\t%s\r\n", line.data());
                    cmdList[i]->execute(argument.data());
                    
                    // Close file, if redirection was requested, then
                    // remove unwanted '\t' (i.e. tab) from that file
                    if(stdOutFile != stdout)
                    {
                        fclose(stdOutFile);
                        CFileUtility::removeTabsFromFile(fileName);
                        cmdList[i]->redirectOutput(stdout);
                    }
                }
            }
        }
        
    } while(!lastLineInFile);

    fprintf(stdOut_, "\r\n");

    fclose(fileHandle);
    return true;
}

bool CCmdLoad::readPythonFile(std::string scriptFileName)
{
    FILE* fileHandle = fopen(scriptFileName.data(), "r");
    
    if(!state_) return false;
    
    if(!fileHandle)
    {
        printf("\r\nError, could not open file: %s!\r\n", scriptFileName.data());
        return false;
    }
    
    fseek(fileHandle, 0L, SEEK_END);
    unsigned long size = ftell(fileHandle);

    fseek(fileHandle, 0L, SEEK_SET);
    std::string fileName;
    fileName.resize(size + 10);
    fread((char*)fileName.data(), sizeof(char), (size_t)size, fileHandle);
    fileName[size] = '\0';
    fclose(fileHandle);

    PyRun_SimpleString(fileName.data());
    return true;
}

bool CCmdLoad::readMassChargeFile(std::string massChargeFileName)
{
    double  Q, m;
    std::string line, text, IDCmp;
    bool lastLineInFile = false;
    FILE* fileHandle = fopen(massChargeFileName.data(), "r");
    int atomIndex;
    
    if(!state_) return false;

    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);

        text = CASCIIUtility::getWord(line, 0);
        
        // If first word is 'ID' then every atom with ID described in
        // the second word will get the charge described in the third word
        if(text == "ID")
        {
            text = CASCIIUtility::getWord(line, 1);
            IDCmp = text;

            text = CASCIIUtility::getWord(line, 2);
            Q = (double)atof(text.data());

            text = CASCIIUtility::getWord(line, 3);
            m = (double)atof(text.data());
            
            for(int i=0; i<(int)state_->atoms_.size(); i++)
            {
                std::string ID = state_->atoms_[i]->getID();
                if(ID == IDCmp)
                {
                    state_->atoms_[i]->Q_ = Q;
                    state_->atoms_[i]->m_ = m;
                }
            }
        }
        
        // If first word is 'AtInd' then every atom with index described in
        // the second word will get the charge described in the third word
        if(text == "AtInd")
        {
            text = CASCIIUtility::getWord(line, 1);
            atomIndex = (int)atoi(text.data());
            
            text = CASCIIUtility::getWord(line, 2);
            Q = (double)atof(text.data());

            text = CASCIIUtility::getWord(line, 3);
            m = (double)atof(text.data());
            
            if(atomIndex < (int)state_->atoms_.size())
            {
                state_->atoms_[atomIndex]->Q_ = Q;
                state_->atoms_[atomIndex]->m_ = m;
            }
        }

    } while(!lastLineInFile);
    
    
    fclose(fileHandle);
    return true;
}

bool CCmdLoad::readQEPosFile(std::string qePosFileName, std::string qeInputFileName)
{
    std::vector<std::string> atomsList;
    std::vector<std::string> speciesList;
    const double auToAngstom = 0.52917721092;
    double X, Y, Z;
    std::string ID, line, text;
    FILE* fileHandlePos = fopen(qePosFileName.data(), "r");
    FILE* fileHandleInput = fopen(qeInputFileName.data(), "r");
    bool lastLineInFile, firstFrameAdded=false;
    int numCoordinates;
    
    
    if(!state_) return false;
    
    if(!fileHandlePos || !fileHandleInput)
    {
        printf("\r\nCould not open files: %s, %s!\r\n", qePosFileName.data(), qeInputFileName.data());
        return false;
    }
    
    
    // Get confiramation from user and purge everything if 'Yes'
    printf("Are you sure you want to delete everything (Yes/No)?: ");
    char stringAnswer[20];
    fgets(stringAnswer, 19, stdin);
    std::string stringYesNo = CASCIIUtility::removeCRLF(stringAnswer);
    if(stringYesNo != "Yes") return true;
    state_->purgeAtomsList();
    

    // Read atomic species first. Atoms in the ATOMIC_POSITION section
    // are labeled according to atoms in the species section ane are
    // in a random order. However, in the output (e.g. *.pos files)
    // the coordinates are listed in order of atomic species.
    do
    {
        line = CFileUtility::readLine(fileHandleInput, lastLineInFile);
        
        if(line.find("ATOMIC_SPECIES") != std::string::npos)
        {
            speciesList.clear();
            
            do
            {
                line = CFileUtility::readLine(fileHandleInput, lastLineInFile);
                if(line.find("/") != std::string::npos) lastLineInFile = true;
                else if(line.find("ATOMIC_POSITIONS") != std::string::npos) lastLineInFile = true;
                else
                {
                    if(!CASCIIUtility::isLineEmpty(line))
                    {
                        ID = CASCIIUtility::getWord(line, 0);
                        speciesList.emplace_back(ID);
                    }
                }
                
            } while(!lastLineInFile);
        }
        
    } while(!lastLineInFile);

    
    // Read, in order, each coordinate species from the ATOMIC_POSITION section
    // and make a corresponding list of all atoms in the correct order (i.e.
    // in the order of the species listed in the ATOMIC_SPECIES section)
    atomsList.clear();
    for(int i=0; i<(int)speciesList.size(); i++)
    {
        fseek(fileHandleInput, 0, SEEK_SET);
        do
        {
            line = CFileUtility::readLine(fileHandleInput, lastLineInFile);
            
            if(line.find("ATOMIC_POSITIONS") != std::string::npos)
            {
                do
                {
                    line = CFileUtility::readLine(fileHandleInput, lastLineInFile);
                    if(line.find("/") != std::string::npos) lastLineInFile = true;
                    else
                    {
                        if(!CASCIIUtility::isLineEmpty(line))
                        {
                            ID = CASCIIUtility::getWord(line, 0);
                            if(ID == speciesList[i])
                                atomsList.emplace_back(ID);
                        }
                    }
                    
                } while(!lastLineInFile);
            }
            
        } while(!lastLineInFile);
    }
    
    
    // Read coordinates *.pos from file
    do
    {
        // Read pos file record header
        line = CFileUtility::readLine(fileHandlePos, lastLineInFile);
        CASCIIUtility::removeWhiteSpace(line);
        numCoordinates = (int)atomsList.size();
        if(!lastLineInFile && (numCoordinates > 0))
        {
            int indexFrame=0;
            
            if(firstFrameAdded) indexFrame = state_->addFrame();
            
            // Read positions
            for(int i=0; i<numCoordinates; i++)
            {
                line = CFileUtility::readLine(fileHandlePos, lastLineInFile);

                text = CASCIIUtility::getWord(line, 0);
                X = atof(text.data()) * auToAngstom;
                
                text = CASCIIUtility::getWord(line, 1);
                Y = atof(text.data()) * auToAngstom;
                
                text = CASCIIUtility::getWord(line, 2);
                Z = atof(text.data()) * auToAngstom;
                
                if(!firstFrameAdded)
                    state_->addAtom(X, Y, Z, atomsList[i].data());
                else
                {
                    if(indexFrame >= 0) state_->setAtomCoordinates(indexFrame, i, X, Y, Z);
                }
                if(i < (int)state_->atoms_.size())
                    state_->atoms_[i]->sigma_ = state_->defaultAtProp_.getWDWRadius(atomsList[i].data());
            }
            
            firstFrameAdded = true;
        }
        
    } while(!lastLineInFile);
    
    fclose(fileHandlePos);
    fclose(fileHandleInput);
    return true;
}
