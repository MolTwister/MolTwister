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

#pragma once
#include "MolTwisterState.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CToolsBase
{
public:
    CToolsBase() = delete;
    CToolsBase(CMolTwisterState* state, FILE* stdOut);

protected:
    CMolTwisterState* state_;
    FILE* stdOut_;
};

END_CUDA_COMPATIBLE()
