#pragma once

class CMolTwisterRestAPI
{
public:
    CMolTwisterRestAPI() = default;

public:
    void run() const;

private:
    static void* threadRun(void*);
};
