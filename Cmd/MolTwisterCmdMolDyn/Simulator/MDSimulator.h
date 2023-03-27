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

// Note! This is an entry point to GPU compiled code.
// Hence, it is necessary to include this file twice
// in different namespaces, within the same file. This
// means that 'pragma once' should not be used

#include <string>
#include "../../../Utilities/CUDAGeneralizations.h"
#include "../../../Utilities/Serializer.h"
#include "../Config/MolDynConfigStruct.h"

BEGIN_CUDA_COMPATIBLE()

class CMDSimulator
{
public:
    CMDSimulator() {}

public:
    static void run(SMolDynConfigStruct config, FILE* stdOut, CSerializer& stateContent);
    static void run(SMolDynConfigStruct config, FILE* stdOut, void* state);
    static void optimize(SMolDynConfigStruct config, FILE* stdOut, CSerializer& stateContent);
    static void optimize(SMolDynConfigStruct config, FILE* stdOut, void* state);
};

END_CUDA_COMPATIBLE()
