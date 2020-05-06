#pragma once
#include <vector>
#include "3DView/MolTwister3DView.h"
#include "Cmd/MolTwisterCmd.h"

class CMolTwister
{
public:
    enum EStartDir { dirWorking=0, dirHome=1, dirInitFile=2 };
    
public:
    CMolTwister();
    ~CMolTwister();
    
public:
    void run(C3DView* view3D);
    std::vector<std::shared_ptr<CAtom>>* getAtomsPtr() { return &state_.atoms_; }
    std::vector<std::shared_ptr<CGLObject>>* getGLObjsPtr() { return &state_.glObjects_; }
    int* getCurrFramePtr() { return &state_.currentFrame_; }
    CDefaultAtomicProperties* getDefaultAtProp() { return &state_.defaultAtProp_; }
    CMolTwisterState* getCurrentState() { return &state_; }
    
protected:
    
private:
    bool _run();
    static void* threadRun(void* arg);
    void initReadline() const;
    static char** commandCompletionFunc(const char* text, int start, int end);
    static char* duplicateString(const char* str);
    static char* commandCompletionGenerator(const char* text, int state);
    std::string readLine() const;
    void readShortcutDirs();
    
public:
    EStartDir startDir_;
    
private:
    CMolTwisterState state_;
    static std::vector<std::shared_ptr<CCmd>> cmdList_;
};
