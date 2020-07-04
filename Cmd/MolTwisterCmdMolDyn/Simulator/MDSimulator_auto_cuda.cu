#include "MDSimulator.h"
#include "../../../MolTwisterState.h"
#include "../MDLoop/MDLoop.h"
#include "../MDLoop/Printf.h"

BEGIN_CUDA_COMPATIBLE()

void CMDSimulator::run(int numMDSteps, int outputStride, int nhChainLength, std::string progressOutFileName, FILE* stdOut, CSerializer& stateContent, int ensemble)
{
    CMolTwisterState state;
    state.serialize(stateContent, false);
    run(numMDSteps, outputStride, nhChainLength, progressOutFileName, stdOut, &state, ensemble);
}

void CMDSimulator::run(int numMDSteps, int outputStride, int nhChainLength, std::string progressOutFileName, FILE* stdOut, void* state, int ensemble)
{
    CMDLoop         MDLoop;
    CSimulationBox  SimBox((CMolTwisterState*)state, stdOut);

    SimBox.bNPTEnsemble = (ensemble == ensembleNPT) ? true : false;
    COut::SetOutputFile(fopen(progressOutFileName.data(), "w"));
    SimBox.InitSystem(nhChainLength);
    MDLoop.RunSimulation(SimBox, numMDSteps, outputStride);
    COut::CloseOutputFile();
}

END_CUDA_COMPATIBLE()
