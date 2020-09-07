#include "MDSimulator.h"
#include "../../../MolTwisterState.h"
#include "../MDLoop/MDLoop.h"
#include "../MDLoop/Printf.h"

BEGIN_CUDA_COMPATIBLE()

void CMDSimulator::run(SMolDynConfigStruct config, FILE* stdOut, CSerializer& stateContent)
{
    CMolTwisterState state;
    state.serialize(stateContent, false);
    run(config, stdOut, &state);
}

void CMDSimulator::run(SMolDynConfigStruct config, FILE* stdOut, void* state)
{
    CMDLoop         MDLoop(config.includeXYZFile_, config.outXYZFile_, config.includeDCDFile_, config.outDCDFile_);
    CSimulationBox  SimBox((CMolTwisterState*)state, stdOut, config);

    COut::SetOutputFile(fopen(config.outInfoFile_.data(), "w"));
    MDLoop.RunSimulation(SimBox, config.numberOfTimeSteps_, config.outputStride_);
    COut::CloseOutputFile();
}

END_CUDA_COMPATIBLE()
