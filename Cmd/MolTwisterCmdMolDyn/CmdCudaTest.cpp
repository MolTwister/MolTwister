#include "CmdCudaTest.h"
#include "../../Utilities/ASCIIUtility.h"
#include "CudaTests/CudaTests_cuda.h"

std::string CCmdCudaTest::getCmd()
{
    return "cudatest";
}

std::vector<std::string> CCmdCudaTest::getCmdLineKeywords()
{
    return { "cudatest" };
}

std::vector<std::string> CCmdCudaTest::getCmdHelpLines()
{
    return {
        "cudatest"
    };
}

std::string CCmdCudaTest::getCmdFreetextHelp()
{
    std::string text;

    text+= "\t\r\n";
    text+= "\t\r\n";
    text+= "\t\r\n";
    text+= "\t\r\n";
    text+= "\t\r\n";
    text+= "\t\r\n";

    return text;
}

std::string CCmdCudaTest::execute(std::vector<std::string>)
{
    lastError_ = "";

    // Test CUDA without the Thrust library
    std::vector<int> A(CCudaTest::getArraySize());
    std::vector<int> B(CCudaTest::getArraySize());

    for(int i=0; i<CCudaTest::getArraySize(); i++)
    {
        A[i] = i;
        B[i] = i+1;
    }
    std::vector<int> cpyA = A;

    CCudaTest::addBIntoA(A.data(), B.data());

    fprintf(stdOut_, "\r\n\tCalculating A_i + B_i, where A_i = i and B_i = i + 1, and i is in [0, %i], using CUDA\r\n", CCudaTest::getArraySize());
    for(int i=0; i<CCudaTest::getArraySize(); i++)
    {
        fprintf(stdOut_, "\t%i + %i = %i\r\n", cpyA[i], B[i], A[i]);
    }

    // Test CUDA with the Thrust library
    fprintf(stdOut_, "\r\n");
    CCudaTest::testModifyAtomList(stdOut_);

    return lastError_;
}

