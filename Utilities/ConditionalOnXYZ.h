#pragma once
#include <stdio.h>
#include <string>
#include <vector>

//////////////////////////////////////////////////////
// Conditionals are of the form x>num, where x can
// be exchanged with y or z. Allowed conditionals are
// <, <=, =, >=, >. Several conditionals can be
// combined using &,* (AND), ^ (XOR), |,+ (OR). Use
// parenthesis to define evaluation order (expressions
// are evaluated from left to right without order,
// other than parenthesis). Use ! prior to any
// conditional or parenthesis to negate (i.e. !=NOT).
// SPACE IS NOT ALLOWED!
//////////////////////////////////////////////////////

class CConditionalOnXYZ
{
private:
    enum ELogic { lAnd=0, lOr, lXOr };
    
private:
    class CConditional
    {
    private:
        enum EConditional { cNone=0, cLess, cLessEq, cEq, cGreaterEq, cGreater, cNotEq };
        enum ECoord { coordX=0, coordY, coordZ };
        
    public:
        CConditional() { condOnXYZ_ = nullptr; conditional_ = cNone; coordinate_ = coordX; condNumber_ = 0.0; isnot_ = false; }
        ~CConditional() { }
        
    public:
        bool parse(std::string conditional, bool isnot);
        bool isFulfilled(double x, double y, double z) const;
        
    private:
        std::shared_ptr<CConditionalOnXYZ> condOnXYZ_;
        EConditional conditional_;
        ECoord coordinate_;
        double condNumber_;
        bool isnot_;
    };

public:
    CConditionalOnXYZ() {}
    CConditionalOnXYZ(std::string conditional) { parse(conditional); }
    ~CConditionalOnXYZ();
    
public:
    bool parse(std::string conditional);
    bool isFulfilled(double x, double y, double z) const;
    static std::string getHelpText(std::string nameOfConditional, bool endWithCRLF);
    
private:
    void purge();
    bool isAND(char c) const { if(c == '&') return true; if(c == '*') return true; return false; }
    bool isXOR(char c) const { if(c == '^') return true; return false; }
    bool isOR(char c) const { if(c == '|') return true; if(c == '+') return true; return false; }
    bool isNOT(char c) const { if(c == '!') return true; return false; }
    std::string getParenthesis(std::string in, int startParenthesis, int& endParenthesisIndex) const;
    
private:
    std::vector<std::shared_ptr<CConditional>> conditionals_;
    std::vector<ELogic> logicBetweenCond_;
};
