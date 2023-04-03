//
// Copyright (C) 2023 Richard Olsen.
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

#include "CmdFFT.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/FFT1D.h"
#include "../../Utilities/FileUtility.h"

std::string CCmdFFT::getCmd()
{
    return "fft";
}

std::vector<std::string> CCmdFFT::getCmdLineKeywords()
{
    return { "fft", "fwd", "rev", "zeropad" };
}

std::vector<std::string> CCmdFFT::getCmdHelpLines()
{
    return {
                "fft <ASCII file with 1 or more columns> <index of columns with real numbers> <index of column with imaginary numbers> <direction> [zeropad]"
           };
}

std::string CCmdFFT::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the FFT (Fast Fourier Transform) of a list of complex numbers, taken from an ASCII file\r\n";
    text+= "\tcontaining several coulumns of data (separated using one or more space characters). The columns are \r\n";
    text+= "\tchosen by their zero base indices. One or two columns can be selected, one for the real part of the\r\n";
    text+= "\tcomplex numners and one for the imaginary part. In case of only real numbers, then the imaginary\r\n";
    text+= "\tcolumn index can be set to -1. The <direction> parameter can be either 'fwd' or 'rev', depending on\r\n";
    text+= "\tif a forward or reverse Fourier transform is to be calculated, respectively.\r\n";
    text+= "\r\n";
    text+= "\tBy invoking 'zeropad' the input dataset is zero padded to achieve a size 2^M, where M is an integer.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. Re Im\r\n";
    text+= "\t2. <transformed real number> <transformed imaginary number>\r\n";
    text+= "\t3. <transformed real number> <transformed imaginary number>\r\n";
    text+= "\t            .\r\n";
    text+= "\t            .\r\n";
    text+= "\t            .\r\n";
    text+= "\tN+1. <transformed real number> <transformed imaginary number>\r\n";
    text+= "\twhere N is the number of entries in the FFT transformed dataset.";

    return text;
}

std::string CCmdFFT::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CFFT1D::EDir dir;
    std::string line, text;
    FILE* file;
    bool lastLineInFile, zeroPad = false;
    int reColumnToFFT, imColumnToFFT;


    //////////////////////////////////////////////////////////////////
    // fwd : x_n = \sum_{m=1}^N y_m e^{-2\pi i m n / N }
    // rev : x_n = \sum_{m=1}^N y_m e^{+2\pi i m n / N }
    //////////////////////////////////////////////////////////////////

    // Open ASCII file with more than one column and select column index
    // or two indices if complex input (first Re and second Im)
    text = CASCIIUtility::getArg(arguments, arg++);
    file = fopen(text.data(), "r");
    if(!file)
    {
        lastError_ = std::string("Error: unable to open file ") + text + std::string("!");
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    reColumnToFFT = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    imColumnToFFT = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "fwd")       dir = CFFT1D::dirFwd;
    else if(text == "rev")  dir = CFFT1D::dirRev;
    else
    {
        lastError_ = std::string("Error: direction should be 'fwd' or 'rev', do not recognize ") + text + std::string("!");
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "zeropad")   zeroPad = true;
    else                    arg--;


    // Read column
    std::vector<CFFT1D::CCplx> aIn;
    do
    {
        CFFT1D::CCplx cplx;

        line = CFileUtility::readLine(file, lastLineInFile);
        text = CASCIIUtility::getWord(line, reColumnToFFT);
        cplx.re_ = atof(text.data());

        if(imColumnToFFT != -1)
        {
            text = CASCIIUtility::getWord(line, imColumnToFFT);
            cplx.im_ = atof(text.data());
        }

        aIn.emplace_back(cplx);

    } while(!lastLineInFile);


    // If requested, zero pad to achieve a dataset size of 2^n to
    // avoid discreete fourier transforms (i.e. avoid not using FFT)
    if(zeroPad) CFFT1D::zeroPad(aIn);
    printf("\r\n\tDataset size (%s) is %i!\r\n", zeroPad ? "zero padded" : "not zero padded", (int)aIn.size());


    // Execute FFT
    CFFT1D fft;
    std::shared_ptr<std::vector<CFFT1D::CCplx>> fftOut = fft.fft1D(aIn, dir);



    // Display result
    printf("\r\n");
    fprintf(stdOut_, "\t%-15s%-15s\r\n", "Re", "Im");
    for(int i=0; i<fftOut->size(); i++)
    {
        fprintf(stdOut_, "\t% -15.8f% -15.8f\r\n", (*fftOut)[i].re_, (*fftOut)[i].im_);
    }


    // Clean up
    fclose(file);

    return lastError_;
}
