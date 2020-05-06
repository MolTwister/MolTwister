#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ConditionalOnXYZ.h"
#include "ASCIIUtility.h"

bool CConditionalOnXYZ::CConditional::parse(std::string conditional, bool isnot)
{
    CASCIIUtility::removeWhiteSpace(conditional);
    
    if(conditional.size() < 1) return false;
    
    if(conditional[0] == '(')
    {
        if(conditional[conditional.size() - 1] != ')') return false;
        
        conditional.erase(conditional.begin());
        if(conditional.size() > 0)
            conditional.erase(conditional.size() - 1);
        
        condOnXYZ_ = std::make_shared<CConditionalOnXYZ>(conditional);
    }
    
    else
    {
        int startNum = 2;
        
        if(conditional.size() < 3) return false;
        
        switch(conditional[0])
        {
            case 'x':
                coordinate_ = coordX;
                break;

            case 'y':
                coordinate_ = coordY;
                break;

            case 'z':
                coordinate_ = coordZ;
                break;
                
            default:
                return false;
        }
        
        if(conditional[1] == '<')
        {
            if(conditional[2] == '=')
            {
                startNum = 3;
                conditional_ = cLessEq;
            }
            else
            {
                conditional_ = cLess;
            }
        }
        else if(conditional[1] == '=')
        {
            conditional_ = cEq;
        }
        else if(conditional[1] == '>')
        {
            if(conditional[2] == '=')
            {
                startNum = 3;
                conditional_ = cGreaterEq;
            }
            else
            {
                conditional_ = cGreater;
            }
        }
        else if((conditional[1] == '!') && (conditional[2] == '='))
        {
            startNum = 3;
            conditional_ = cNotEq;
        }
        else
        {
            return false;
        }
        
        std::string condNumber;
        if(startNum >= conditional.size()) return false;
        for(int i=startNum; i<conditional.size(); i++)
            condNumber+= conditional[i];
        
        condNumber_ = atof(condNumber.data());
    }
    
    isnot_ = isnot;
    
    return true;
}

bool CConditionalOnXYZ::CConditional::isFulfilled(double x, double y, double z) const
{
    bool result;
    
    if(condOnXYZ_)
    {
        result = condOnXYZ_->isFulfilled(x, y, z);
    }
    else if(conditional_ != cNone)
    {
        double coord;
        
        if(coordinate_ == coordX) coord = x;
        else if(coordinate_ == coordY) coord = y;
        else if(coordinate_ == coordZ) coord = z;
        else return false;
        
        if(conditional_ == cLess)              result = (coord < condNumber_) ? true : false;
        else if(conditional_ == cLessEq)       result = (coord <= condNumber_) ? true : false;
        else if(conditional_ == cEq)           result = (fabs(coord - condNumber_) < double(FLT_MIN)) ? true : false;
        else if(conditional_ == cGreater)      result = (coord > condNumber_) ? true : false;
        else if(conditional_ == cGreaterEq)    result = (coord >= condNumber_) ? true : false;
        else if(conditional_ == cNotEq)        result = (fabs(coord - condNumber_) > double(FLT_MIN)) ? true : false;
        else return false;
    }
    else
    {
        return false;
    }
    
    if(isnot_) result = result ? false : true;
    
    return result;
}



CConditionalOnXYZ::~CConditionalOnXYZ()
{
    purge();
}

void CConditionalOnXYZ::purge()
{
    conditionals_.clear();
    logicBetweenCond_.clear();
}

std::string CConditionalOnXYZ::getParenthesis(std::string in, int startParenthesis, int& endParenthesisIndex) const
{
    int parNum = 0;
    std::string parenthesis;
       
    for(int i=startParenthesis; i<in.size(); i++)
    {
        if(in[i] == '(') parNum++;
        if(in[i] == ')') parNum--;
        
        parenthesis+= in[i];
        
        if(parNum <= 0)
        {
            endParenthesisIndex = i;
            return parenthesis;
        }
    }

    endParenthesisIndex = (int)in.size() - 1;
    return parenthesis;
}

bool CConditionalOnXYZ::parse(std::string conditional)
{
    std::string text;
    bool isnot = false;
    
    purge();
    for(int i=0; i<conditional.size(); i++)
    {
        if(isNOT(conditional[i]))
        {
            isnot = true;
        }
        else if(isAND(conditional[i]))
        {
            auto cond = std::make_shared<CConditional>();
            cond->parse(text, isnot);
            text.clear();
            isnot = false;
            conditionals_.emplace_back(cond);
            logicBetweenCond_.emplace_back(lAnd);
        }
        else if(isOR(conditional[i]))
        {
            auto cond = std::make_shared<CConditional>();
            cond->parse(text, isnot);
            text.clear();
            isnot = false;
            conditionals_.emplace_back(cond);
            logicBetweenCond_.emplace_back(lOr);
        }
        else if(isXOR(conditional[i]))
        {
            auto cond = std::make_shared<CConditional>();
            cond->parse(text, isnot);
            text.clear();
            isnot = false;
            conditionals_.emplace_back(cond);
            logicBetweenCond_.emplace_back(lXOr);
        }
        else if(conditional[i] == '(')
        {
            text = getParenthesis(conditional, i, i);
        }
        else
        {
            text+= conditional[i];
        }

        if(i == (conditional.size() - 1))
        {
            auto cond = std::make_shared<CConditional>();
            cond->parse(text, isnot);
            text.clear();
            isnot = false;
            conditionals_.emplace_back(cond);
        }
    }
    
    return true;
}

bool CConditionalOnXYZ::isFulfilled(double x, double y, double z) const
{
    bool ret = false;
    
    for(int i=0; i<conditionals_.size(); i++)
    {
        if(i == 0)
        {
            ret = conditionals_[i]->isFulfilled(x, y, z);
        }
        else
        {
            if((i-1) >= logicBetweenCond_.size()) return false;
            if(logicBetweenCond_[i-1] == lAnd)
            {
                ret&= conditionals_[i]->isFulfilled(x, y, z);
            }
            else if(logicBetweenCond_[i-1] == lOr)
            {
                ret|= conditionals_[i]->isFulfilled(x, y, z);
            }
            else if(logicBetweenCond_[i-1] == lXOr)
            {
                ret^= conditionals_[i]->isFulfilled(x, y, z);
            }
            else
            {
                return false;
            }
        }
    }
    
    return ret;
}

std::string CConditionalOnXYZ::getHelpText(std::string nameOfConditional, bool endWithCRLF)
{
    std::string text;

    text+= std::string("\tThe ") + nameOfConditional + std::string(" parameter consists of a string without space, where several\r\n");
    text+= "\tconditions are specified, separated by boolean operators. The boolean\r\n";
    text+= "\toperators that can be applied are:\r\n";
    text+= "\t* not: !\r\n";
    text+= "\t* and: &\r\n";
    text+= "\t* and: *\r\n";
    text+= "\t* or: |\r\n";
    text+= "\t* or: +\r\n";
    text+= "\t* xor: ^\r\n";
    text+= "\tEach condition is on x, y or z and can be of the following types\r\n";
    text+= "\t* equal: =\r\n";
    text+= "\t* not equal: !=\r\n";
    text+= "\t* greater than: >\r\n";
    text+= "\t* greater than or equal: >=\r\n";
    text+= "\t* less than: <\r\n";
    text+= "\t* less than or equal: <=\r\n";
    text+= "\tFor example,\r\n";
    text+= "\t* x>5\r\n";
    text+= "\t* x>5&x<10\r\n";
    text+= "\t* (x>5&x<10)|(x>20&x<23)\r\n";
    text+= "\t* !(x>5&x<10)";

    if(endWithCRLF) text+= "\r\n";
    return text;
}
