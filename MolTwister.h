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

#pragma once
#include <vector>
#include "MolTwisterTutorialPool.h"
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
    CMolTwisterTutorialPool tutorialPool_;
};
