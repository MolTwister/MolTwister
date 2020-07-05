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
    CMDLoop         MDLoop;
    CSimulationBox  SimBox((CMolTwisterState*)state, stdOut);

    SimBox.bNPTEnsemble = (config.ensemble_ == SMolDynConfigStruct::ensembleNPT) ? true : false;
    COut::SetOutputFile(fopen(config.outInfoFile_.data(), "w"));
    SimBox.InitSystem(config.temperatureNHChainLength_, config.pressureNHChainLength_);
    MDLoop.RunSimulation(SimBox, config.numberOfTimeSteps_, config.outputStride_);
    COut::CloseOutputFile();
}

END_CUDA_COMPATIBLE()
