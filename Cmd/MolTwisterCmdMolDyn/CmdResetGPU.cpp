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

