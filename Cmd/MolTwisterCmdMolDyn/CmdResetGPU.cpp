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

#include "CmdResetGPU.h"
#include "../../Utilities/ASCIIUtility.h"
#include "Cmd/Tools/CudaDeviceList.h"

std::string CCmdResetGPU::getCmd()
{
    return "resetgpu";
}

std::vector<std::string> CCmdResetGPU::getCmdLineKeywords()
{
    return { "resetgpu" };
}

std::vector<std::string> CCmdResetGPU::getCmdHelpLines()
{
    return {
        "resetgpu <gpu index (zero based)>"
    };
}

std::string CCmdResetGPU::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tWill perform a soft reset of the selected gpu. To list the available\r\n";
    text+= "\tGPUs on the system, use the 'get gpuinfo' command. This will list the\r\n";
    text+= "\tGPUs in order of their assigned indices, with index 0 being the firs.";

    return text;
}

std::string CCmdResetGPU::execute(std::vector<std::string> arguments)
{
    lastError_ = "";
    size_t arg = 0;

    int deviceToReset = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    CCudaDeviceList deviceList;
    int activeDevice = deviceList.GetActiveDevice();
    deviceList.SetActiveDevice(deviceToReset);
    deviceList.ResetActiveGPU();
    deviceList.SetActiveDevice(activeDevice);

    return lastError_;
}

