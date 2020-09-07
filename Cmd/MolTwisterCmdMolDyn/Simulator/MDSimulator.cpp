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
    CMDLoop mdLoop(config.includeXYZFile_, config.outXYZFile_, config.includeDCDFile_, config.outDCDFile_);
    CSimulationBox simBox((CMolTwisterState*)state, stdOut, config);

    COut::setOutputFile(fopen(config.outInfoFile_.data(), "w"));
    mdLoop.runSimulation(simBox, config.numberOfTimeSteps_, config.outputStride_);
    COut::closeOutputFile();
}

END_CUDA_COMPATIBLE()
