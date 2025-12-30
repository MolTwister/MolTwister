//
// Copyright (C) 2025 Richard Olsen.
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

#pragma once
#include <stdio.h>
#include <string>
#include <vector>
#include <tuple>
#include <functional>
#include "3DVector.h"
#include "3DRect.h"

BEGIN_CUDA_COMPATIBLE()

class CLammpsTrjFile
{
public:
    class CAtom
    {
    public:
        CAtom() = default;
        CAtom(const std::string& line, const int& indexAtom, const int& indexType, const int& indexX, const int& indexY, const int& indexZ);

    public:
        void unscaleX(const std::pair<double, double>& bbXSpan);
        void unscaleY(const std::pair<double, double>& bbYSpan);
        void unscaleZ(const std::pair<double, double>& bbZSpan);
        int getAtomIndex() const;
        C3DVector getPosition() const;

    private:
        bool parse(const std::string& line, const int& indexAtom, const int& indexType, const int& indexX, const int& indexY, const int& indexZ);

    private:
        int atomNumber_ = -1;
        int typeIndex_ = -1;
        double x_ = 0.0;
        double y_ = 0.0;
        double z_ = 0.0;
    };

    class CRecord
    {
    public:
        CRecord() = default;
        CRecord(std::shared_ptr<std::ifstream> file);

    public:
        bool read(std::shared_ptr<std::ifstream> file);
        int getNumAtoms() const;
        int getTimeStepIndex() const;

        C3DVector getCoordinate(const int& index) const;
        C3DRect getPBC() const;

    private:
        int timeStep_ = -1;
        int numAtoms_ = -1;
        std::pair<double, double> pbcXSpan_;
        std::pair<double, double> pbcYSpan_;
        std::pair<double, double> pbcZSpan_;
        std::vector<CAtom> atoms_;
    };

public:
    CLammpsTrjFile() = default;
    ~CLammpsTrjFile();

public:
    bool open(std::string fileName);
    void close();
    void gotoRecord(int recordIndex);
    int getNumRecords() const;
    int getNumCoordinatesInRecord() const;
    int getTimeStepIndex() const;
    C3DVector getCoordinate(int atomIndex) const;
    C3DRect getCurrentPBC() const;

private:
    std::shared_ptr<std::ifstream> file_;
    std::vector<std::streampos> recordFilePositions_;
    std::shared_ptr<CRecord> currentRecord_ = nullptr;
};

END_CUDA_COMPATIBLE()
