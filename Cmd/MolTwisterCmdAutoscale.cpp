#include <iostream>
#include "MolTwisterCmdAutoscale.h"

void CCmdAutoscale::onAddKeywords()
{
    addKeyword("autoscale");
}

std::string CCmdAutoscale::getHelpString() const
{ 
    std::string text;
    
    text+= "\tUsage: autoscale\r\n";
    text+= "\r\n";
    text+= "\tAutoscale the 3D view by searching for the boundary of all atoms inside the\r\n";
    text+= "\t3D view and subsequently adjusting the viewing frustum to these boundaries.\r\n";
    text+= "\tThe camera is placed along the x-axis and will point towards the center of\r\n";
    text+= "\tthe located atomic boundaries.\r\n";
    
    return text;
}

void CCmdAutoscale::execute(std::string)
{
    if(!state_) return;

    if(state_->view3D_) state_->view3D_->requestUpdate(true);
}
