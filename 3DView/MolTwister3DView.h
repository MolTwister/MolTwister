#pragma once
#include "Utilities/3DVector.h"
#include "Utilities/3DRect.h"
#include "Utilities/CUDAGeneralizations.h"
#include "Utilities/Serializer.h"
#include "DefaultAtomicProperties.h"
#include "ExpLookup.h"
#include "MolTwisterAtom.h"
#include "MolTwisterGLObject.h"
#include <memory>

BEGIN_CUDA_COMPATIBLE()

class C3DView
{
private:
    enum ESelModes { selmodeNone=0, selmodeAtoms=1 };
    
private:
    class CCamera
    {
    public:
        CCamera() { pos_ = C3DVector(3.0, 0.0, 0.0); lookAt_ = C3DVector(0.0, 0.0, 0.0); up_ = C3DVector(0.0, 0.0, 1.0); zoomFactor_ = 1.0; }

    public:
        void serialize(CSerializer& io, bool saveToStream);

    public:
        C3DVector pos_;
        C3DVector lookAt_;
        C3DVector up_;
        double zoomFactor_;
    };
    
    class CScreenCoord
    {
    public:
        CScreenCoord() { x_ = y_ = 0; }
        CScreenCoord(int X, int Y) { x_ = X; y_ = Y; }
        CScreenCoord(const CScreenCoord& src) { x_ = src.x_; y_ = src.y_; }

    public:
        void serialize(CSerializer& io, bool saveToStream);
        double length2() const { return double(x_*x_ + y_*y_); }
        CScreenCoord operator-(const CScreenCoord& rhs) const { CScreenCoord ret(x_ - rhs.x_, y_ - rhs.y_); return ret; }
        CScreenCoord operator+(const CScreenCoord& rhs) const { CScreenCoord ret(x_ + rhs.x_, y_ + rhs.y_); return ret; }
        CScreenCoord operator=(const CScreenCoord& src) { x_ = src.x_; y_ = src.y_; return *this; }
        
    public:
        int x_;
        int y_;
    };
    
    class CResolution
    {
    public:
        CResolution() { sphere_ = 24; cylinder_ = 36; }

    public:
        void serialize(CSerializer& io, bool saveToStream);

    public:
        int sphere_;
        int cylinder_;
    };
    
    class CArg
    {
    public:
        CArg() { argc_ = 0; argv_ = nullptr; }
        CArg(int argc, char *argv[]) { argc_ = argc; argv_ = argv; deleteArgs_ = false; }
        ~CArg();

    public:
        void serialize(CSerializer& io, bool saveToStream);

    public:
        int argc_;
        char** argv_;

    private:
        bool deleteArgs_;
    };
    
public:
    C3DView(int argc, char *argv[]);
    
public:
    void serialize(CSerializer& io, bool saveToStream,
                   std::vector<std::shared_ptr<CAtom>>* atoms=nullptr, std::vector<std::shared_ptr<CGLObject>>* glObjects=nullptr,
                   int* currentFrame=nullptr, CDefaultAtomicProperties* defAtProp=nullptr);

    void show(std::vector<std::shared_ptr<CAtom>>* atoms, std::vector<std::shared_ptr<CGLObject>>* glObjects, int* currentFrame, CDefaultAtomicProperties* defAtProp);
    void requestUpdate(bool updateCameraPos) { if(updateCameraPos) updateRequested_ = 2; else updateRequested_ = 1; }
    void requestFullScreen(bool on) { if(on) fullscreenRequested_ = 1; else fullscreenRequested_ = 2; }
    void setOrthoView(bool set=true) { orthoView_ = set; }
    static void setRedrawLimit(int numAtoms) { numAtomsBeforeNoDraw_ = numAtoms; }
    void setFog(bool enabled) { fogEnabled_ = enabled; }
    void enableViewOfBondsAcrossPBC(bool set=true) { viewBondsAcrossPBC_ = set; }
    void requestQuit() { requestQuit_ = true; }
    C3DRect getPBC() const { return pbc_; }
    static C3DRect calcPBC(int frame=-1);
    void enableUserPBC(bool enable=true) { applyUserDefPBC_ = enable; }
    bool isUserPBCEnabled() const { return applyUserDefPBC_; }
    void setUserPBC(C3DRect pbc) { pbcUser_ = pbc; }
    C3DRect getUserPBC() const { return pbcUser_; }
    
protected:
    static void onRender();
    static void onReshape(int w, int h);
    static void onTimer(int value);
    static void onMouseClick(int button, int state, int x, int y);
    static void onMouseMove(int x, int y);
    static void onKeyboard(unsigned char key, int x, int y);
    static void onSpecialFunc(int key, int x, int y);
    
private:
    static void initOpenGL();
    static void initScene();
    static void initSelColorRot();
    static void drawAtom(C3DVector R, double radius, float r, float g, float b, bool selection=false, C3DVector selColor=C3DVector(1.0, 1.0, 1.0));
    static bool drawBond(C3DVector R1, C3DVector R2, double radius, float r, float g, float b);
    static double drawBitmapText(const char* text, void* glutBitmapFont, double x, double y, double r, double g, double b);
    static double drawBitmapText(const char* text, void* glutBitmapFont, double x, double y, C3DVector color);
    static void drawPBCGrid(ESelModes mode);
    static void drawMolecules(ESelModes mode);
    static void drawGLObjects(ESelModes mode);
    static void drawSelectionHUD(ESelModes mode);
    static void drawIsoSurfaces(ESelModes mode);
    static double calcVDWDensity(double x, double y, double z);
    static void drawScene(ESelModes mode);
    static void update(bool updateCameraPos);
    static void recalcFrustum(double& frustLow, double& frustHigh, ESelModes mode=selmodeNone);
    static void pickAtoms(int x, int y);
    static C3DVector currFrmAtVec(const CAtom& atom);
    
private:
    static int updateRequested_;
    static int fullscreenRequested_;
    static bool requestQuit_;
    static bool orthoView_;
    static bool viewAxes_;
    static bool viewIsoSurface_;
    static bool viewBondsAcrossPBC_;
    static std::vector<std::shared_ptr<CAtom>>* atoms_;
    static std::vector<std::shared_ptr<CGLObject>>* glObjects_;
    static int* currentFrame_;
    static CCamera camera_;
    static bool leftMButtonPressed_;
    static bool middleMButtonPressed_;
    static bool rightMButtonPressed_;
    static CScreenCoord coordLastClick_;
    static CScreenCoord lastWindowSize_;
    static CResolution primitiveRes_;
    static CArg progArg_;
    static std::vector<C3DVector> selColorRot_;
    static CDefaultAtomicProperties* defaultAtProp_;
    static CExpLookup expLookup_;
    static bool applyUserDefPBC_;
    static C3DRect pbcUser_;
    static C3DRect pbc_;
    static int numAtomsBeforeNoDraw_;
    static bool fogEnabled_;
};

END_CUDA_COMPATIBLE()
