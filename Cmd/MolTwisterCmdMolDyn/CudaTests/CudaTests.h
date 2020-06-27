#pragma once
#include <stdio.h>

class CCudaTest
{
public:
    static void addBIntoA(int* A, int* B);
    static void testModifyAtomList(FILE* stdOut);
    static int getArraySize();
};
