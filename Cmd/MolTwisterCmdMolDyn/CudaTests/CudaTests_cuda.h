#pragma once

class CCudaTest_cuda
{
public:
    static void addBIntoA(int* A, int* B);
    static void testModifyAtomList(FILE* stdOut);
    static int getArraySize();
};
