#pragma once

class CBashControl
{
public:
    CBashControl() {}
    
public:
    static void saveCursorPos();
    static void restoreCursorPos();
    static void hideCursor();
    static void showCursor();
    static void clearLastLine();
    static void moveCursorUp();
    static void moveCursorDown();
    static void moveCursorTo(int x, int y);
};

class CBashColor
{
public:
    enum EColor { colNone = 0, colBlack, colRed, colGreen, colYellow, colBlue, colMagenta, colCyan, colWhite };
    enum ESpecial { specReset = 0, specBright, specDim, specUnderscore, specBlink, specReverese, SpecHidden };
    
public:
    CBashColor() {}
    
public:
    static void setColor(EColor foreground = colNone, EColor background = colNone);
    static void setSpecial(ESpecial special = specReset);
    static void clearScreen();
};
