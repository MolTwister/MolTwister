#include <iostream>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>
#include <float.h>
#include "Utilities/ASCIIUtility.h"
#include "MarchingCubes.h"
#include "Cmd/Tools/MolTwisterStateTools.h"
#include "MolTwister3DView.h"

int C3DView::updateRequested_ = 0;
int C3DView::fullscreenRequested_ = 0;
bool C3DView::requestQuit_ = false;
bool C3DView::orthoView_ = false;
bool C3DView::viewAxes_ = true;
bool C3DView::viewIsoSurface_ = false;
bool C3DView::viewBondsAcrossPBC_ = true;
std::vector<std::shared_ptr<CAtom>>* C3DView::atoms_ = nullptr;
std::vector<std::shared_ptr<CGLObject>>* C3DView::glObjects_ = nullptr;
int* C3DView::currentFrame_ = nullptr;
C3DView::CCamera C3DView::camera_ = C3DView::CCamera();
bool C3DView::leftMButtonPressed_ = false;
bool C3DView::middleMButtonPressed_ = false;
bool C3DView::rightMButtonPressed_ = false;
C3DView::CScreenCoord C3DView::coordLastClick_ = C3DView::CScreenCoord();
C3DView::CScreenCoord C3DView::lastWindowSize_ = C3DView::CScreenCoord();
C3DView::CResolution C3DView::primitiveRes_ = C3DView::CResolution();
C3DView::CArg C3DView::progArg_ = C3DView::CArg();
int C3DView::numSelRotColors_ = 0;
C3DVector C3DView::selColorRot_[50] = {};
CDefaultAtomicProperties* C3DView::defaultAtProp_ = nullptr;
CExpLookup C3DView::expLookup_;
bool C3DView::applyUserDefPBC_ = false;
C3DRect C3DView::pbcUser_ = C3DRect();
C3DRect C3DView::pbc_ = C3DRect();
int C3DView::numAtomsBeforeNoDraw_ = 7000;
bool C3DView::fogEnabled_ = false;

C3DView::C3DView(int argc, char *argv[])
{
    progArg_ = CArg(argc, argv);

    updateRequested_ = 0;
    fullscreenRequested_ = 0;
    requestQuit_ = false;
    orthoView_ = false;
    viewAxes_ = true;
    viewIsoSurface_ = false;
    viewBondsAcrossPBC_ = true;
    atoms_ = nullptr;
    glObjects_ = nullptr;
    currentFrame_ = nullptr;
    applyUserDefPBC_ = false;
    pbcUser_.rLow_ = C3DVector(-1.0, -1.0, -1.0);
    pbcUser_.rHigh_ = C3DVector(1.0, 1.0, 1.0);
    pbc_.rLow_ = C3DVector(-1.0, -1.0, -1.0);
    pbc_.rHigh_ = C3DVector(1.0, 1.0, 1.0);
 
    leftMButtonPressed_ = false;
    middleMButtonPressed_ = false;
    rightMButtonPressed_ = false;
    
    primitiveRes_.sphere_ = 10;
    primitiveRes_.cylinder_ = 10;
    
    defaultAtProp_ = nullptr;

    expLookup_.init(-10.0, 10.0, 10000);
    initSelColorRot();
}

C3DVector C3DView::currFrmAtVec(const CAtom &atom)
{
    if(*currentFrame_ >= static_cast<int>(atom.r_.size())) return C3DVector(0.0, 0.0, 0.0);
    
    return atom.r_[static_cast<size_t>(*currentFrame_)];
}

C3DRect C3DView::calcPBC(int frame)
{
    C3DRect pbc;
    
    if(applyUserDefPBC_)
    {
        pbc = pbcUser_;
    }
    
    else
    {
        if(atoms_ && (atoms_->size() > 0))
        {
            const double fltMax = static_cast<double>(FLT_MAX);
            pbc.rLow_ = C3DVector(fltMax, fltMax, fltMax);
            pbc.rHigh_ = C3DVector(-fltMax, -fltMax, -fltMax);
            
            C3DVector r;
            for(size_t i=0; i<atoms_->size(); i++)
            {
                if(currentFrame_ && (frame < 0)) r = currFrmAtVec(*(*atoms_)[i]);
                else
                {
                    if((frame >= 0) && (frame < static_cast<int>((*atoms_)[i]->r_.size()))) r = (*atoms_)[i]->r_[static_cast<size_t>(frame)];
                    else
                    {
                        pbc.rLow_ = C3DVector(-1.0, -1.0, -1.0);
                        pbc.rHigh_ = C3DVector(1.0, 1.0, 1.0);
                        return pbc;
                    }
                }
                
                if(r.x_ < pbc.rLow_.x_) pbc.rLow_.x_ = r.x_;
                if(r.y_ < pbc.rLow_.y_) pbc.rLow_.y_ = r.y_;
                if(r.z_ < pbc.rLow_.z_) pbc.rLow_.z_ = r.z_;
                
                if(r.x_ > pbc.rHigh_.x_) pbc.rHigh_.x_ = r.x_;
                if(r.y_ > pbc.rHigh_.y_) pbc.rHigh_.y_ = r.y_;
                if(r.z_ > pbc.rHigh_.z_) pbc.rHigh_.z_ = r.z_;
            }
        }
        else
        {
            pbc.rLow_ = C3DVector(-1.0, -1.0, -1.0);
            pbc.rHigh_ = C3DVector(1.0, 1.0, 1.0);
        }
    }

    if(pbc.getWidthX() < 1E-5) { pbc.rLow_.x_-= 1.0; pbc.rHigh_.x_+= 1.0; }
    if(pbc.getWidthY() < 1E-5) { pbc.rLow_.y_-= 1.0; pbc.rHigh_.y_+= 1.0; }
    if(pbc.getWidthZ() < 1E-5) { pbc.rLow_.z_-= 1.0; pbc.rHigh_.z_+= 1.0; }

    return pbc;
}

void C3DView::initOpenGL()
{    
    glutInit(&progArg_.argc_, progArg_.argv_);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    int w = glutGet(GLUT_SCREEN_WIDTH);
    int h = glutGet(GLUT_SCREEN_HEIGHT);
    glutInitWindowSize(w/2, h);
    glutInitWindowPosition(w/2, 0);
    glutCreateWindow("3D View");
    
    glutDisplayFunc(onRender);
    glutReshapeFunc(onReshape);
    glutMouseFunc(onMouseClick);
    glutMotionFunc(onMouseMove);
    glutKeyboardFunc(onKeyboard);
    glutSpecialFunc(onSpecialFunc);
    glutTimerFunc(100, onTimer, 0);
}

void C3DView::initScene()
{
    float light0Pos[] = { 100.0f, 100.0f, 100.0f, 0.0f };
    float light0Diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    float light0Specular[] = { 1.0f, 1.0f, 1.0f, 0.0f };
    float lightAmbient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    float matSpecular[] = { 1.0f, 1.0f, 1.0f, 0.15f };
    float matShininess[] = { 50.0f };
    
    // Generic setup + clear color
    float backColor[] = { 0.3f, 0.3f, 0.3f };
    glEnable(GL_DEPTH_TEST);
    glClearColor(backColor[0], backColor[1], backColor[2], 0.0f);
    glShadeModel(GL_SMOOTH);
    glDepthRange(0.0, 1.0);
    
    // Setup lighting
    glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0Diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0Specular);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lightAmbient);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Setup fog
    float fogColor[] = { backColor[0], backColor[1], backColor[2], 1.0f };
    int fogMode = GL_LINEAR;
    glFogi(GL_FOG_MODE, fogMode);
    glFogfv(GL_FOG_COLOR, fogColor);
    glHint(GL_FOG_HINT, GL_DONT_CARE);

    // Define material to be used for front and back faces
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matSpecular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, matShininess);
}

void C3DView::initSelColorRot()
{
    int index = 0;
    
    selColorRot_[index++] = C3DVector(0.4, 0.4, 0.6);
    selColorRot_[index++] = C3DVector(0.4, 0.6, 0.4);
    selColorRot_[index++] = C3DVector(0.4, 0.6, 0.6);
    selColorRot_[index++] = C3DVector(0.6, 0.4, 0.4);
    selColorRot_[index++] = C3DVector(0.6, 0.4, 0.6);
    selColorRot_[index++] = C3DVector(0.6, 0.6, 0.4);
    selColorRot_[index++] = C3DVector(0.6, 0.6, 0.6);

    selColorRot_[index++] = C3DVector(0.0, 0.0, 0.6);
    selColorRot_[index++] = C3DVector(0.0, 0.6, 0.0);
    selColorRot_[index++] = C3DVector(0.0, 0.6, 0.6);
    selColorRot_[index++] = C3DVector(0.6, 0.0, 0.0);
    selColorRot_[index++] = C3DVector(0.6, 0.0, 0.6);
    selColorRot_[index++] = C3DVector(0.6, 0.6, 0.0);
    selColorRot_[index++] = C3DVector(0.6, 0.6, 0.6);

    selColorRot_[index++] = C3DVector(0.0, 0.0, 1.0);
    selColorRot_[index++] = C3DVector(0.0, 1.0, 0.0);
    selColorRot_[index++] = C3DVector(0.0, 1.0, 1.0);
    selColorRot_[index++] = C3DVector(1.0, 0.0, 0.0);
    selColorRot_[index++] = C3DVector(1.0, 0.0, 1.0);
    selColorRot_[index++] = C3DVector(1.0, 1.0, 0.0);
    selColorRot_[index++] = C3DVector(1.0, 1.0, 1.0);

    selColorRot_[index++] = C3DVector(0.4, 0.4, 1.0);
    selColorRot_[index++] = C3DVector(0.4, 1.0, 0.4);
    selColorRot_[index++] = C3DVector(0.4, 1.0, 1.0);
    selColorRot_[index++] = C3DVector(1.0, 0.4, 0.4);
    selColorRot_[index++] = C3DVector(1.0, 0.4, 1.0);
    selColorRot_[index++] = C3DVector(1.0, 1.0, 0.4);
    selColorRot_[index++] = C3DVector(1.0, 1.0, 1.0);
    
    numSelRotColors_ = index;
}

void C3DView::drawAtom(C3DVector R, double radius, float r, float g, float b, bool selection, C3DVector selColor)
{
    float color[] = { r, g, b, 1.0 };

    glTranslatef((float)R.x_, (float)R.y_, (float)R.z_);

    if(selection)
    {
        float colorSel[] = { (float)selColor.x_, (float)selColor.y_, (float)selColor.z_, 0.6f };

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorSel);
        glEnable(GL_BLEND);
        glDepthMask(GL_FALSE);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        glutSolidSphere(radius * 1.3, primitiveRes_.sphere_, primitiveRes_.sphere_);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
    }
    else
    {
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
        glutSolidSphere(radius, primitiveRes_.sphere_, primitiveRes_.sphere_);
    }

    glTranslatef(-(float)R.x_, -(float)R.y_, -(float)R.z_);
}

bool C3DView::drawBond(C3DVector R1, C3DVector R2, double radius, float r, float g, float b)
{
    const double twoPi = 2.0*M_PI;
    const double deltaTheta = twoPi / double(primitiveRes_.cylinder_);
    C3DVector L, u, v, na, nb;
    C3DVector p1a, p2a;
    C3DVector p1b, p2b;
    float color[] = { r, g, b, 1.0 };
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);

    // Calculate unit vector u, used as u-axis for base disc of cylinder
    L = R2 - R1;
    if(L.z_ != 0.0)
    {
        u.x_ = u.y_ = 1.0;
        u.z_ = -(L.x_ + L.y_) / L.z_;
    }
    else if(L.y_ != 0.0)
    {
        u.x_ = u.z_ = 1.0;
        u.y_ = -(L.x_ + L.z_) / L.y_;
    }
    else if(L.x_ != 0.0)
    {
        u.y_ = u.z_ = 1.0;
        u.x_ = -(L.y_ + L.z_) / L.x_;
    }
    else
    {
        return false;
    }
    u.normalize();
    
    // Calculate unit vector v, used as v-axis for base disc of cylinder
    v = L.cross(u);
    v.normalize();
    
    // Calculate lines p1(theta) to p2(theta) along cylinder at angle theta
    for(double theta=0; theta<=twoPi; theta+=deltaTheta)
    {
        na = u*cos(theta) + v*sin(theta);
        p1a = R1 + na*radius;
        p2a = p1a + L;
        
        if(theta > 1E-4)
        {
            glBegin(GL_POLYGON);
            glNormal3f((float)na.x_, (float)na.y_, (float)na.z_);
            glVertex3f((float)p1a.x_, (float)p1a.y_, (float)p1a.z_);
            glNormal3f((float)na.x_, (float)na.y_, (float)na.z_);
            glVertex3f((float)p2a.x_, (float)p2a.y_, (float)p2a.z_);
            glNormal3f((float)nb.x_, (float)nb.y_, (float)nb.z_);
            glVertex3f((float)p2b.x_, (float)p2b.y_, (float)p2b.z_);
            glNormal3f((float)nb.x_, (float)nb.y_, (float)nb.z_);
            glVertex3f((float)p1b.x_, (float)p1b.y_, (float)p1b.z_);
            glEnd();
        }

        nb = na;
        p1b = p1a;
        p2b = p2a;
    }
    
    return true;
}

double C3DView::drawBitmapText(const char* text, void* glutBitmapFont, double x, double y, double r, double g, double b)
{
    const char* p = text;
    double lastXRasterPos=0.0;
    int rasterPos[4];
    int viewport[4];
    
    
    // Save current projection matrix and clear the new
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    
    // Save the current model view matrix and clear the new
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    // Disable depth testing to ensure top rendering is possible,
    // set color (only possible if no lighting effects, and must
    // be done before glRasterpos, where glColor sets raster color,
    // else glColor will set object color)
    glDisable(GL_DEPTH_TEST); 
    
    glDisable(GL_LIGHTING);
    glColor4f((float)r, (float)g, (float)b, 1.0f);
    glRasterPos2f((float)x, (float)y);
    
    // Draw characters of text
    do
    {
        glutBitmapCharacter(glutBitmapFont, *p);
        
    } while( *(++p));
    
    glGetIntegerv(GL_CURRENT_RASTER_POSITION, rasterPos);
    glGetIntegerv(GL_VIEWPORT, viewport);
    if(viewport[2] != 0)
    {
        lastXRasterPos = 2.0 * double(rasterPos[0]) / double(viewport[2]) - 1.0;
    }
    
    // Restore previous light, depth, projection and view states
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    
    return lastXRasterPos;
}

double C3DView::drawBitmapText(const char* text, void* glutBitmapFont, double x, double y, C3DVector color)
{
    return drawBitmapText(text, glutBitmapFont, x, y, color.x_, color.y_, color.z_);
}

void C3DView::drawPBCGrid(ESelModes mode)
{
    if(!viewAxes_) return;
    
    double  stepX = pbc_.getWidthX() / 10.0;
    double  stepY = pbc_.getWidthY() / 10.0;
    double  stepZ = pbc_.getWidthZ() / 10.0;
    double  maxSize = stepX;
    double  textGap;

    if(stepY > maxSize) maxSize = stepY;
    if(stepZ > maxSize) maxSize = stepZ;
    
    if(mode == selmodeAtoms) return;

    // Draw grid planes
    glDisable(GL_LIGHTING);
    glColor4f(0.4f, 0.4f, 0.4f, 0.0f);
    
    for(double x=pbc_.rLow_.x_; x<=(pbc_.rHigh_.x_+stepX/2.0); x+=stepX)
    {
        glBegin(GL_LINES);
        glVertex3f((float)x, (float)pbc_.rLow_.y_, (float)pbc_.rLow_.z_);
        glVertex3f((float)x, (float)pbc_.rLow_.y_, (float)pbc_.rHigh_.z_);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f((float)x, (float)pbc_.rLow_.y_, (float)pbc_.rLow_.z_);
        glVertex3f((float)x, (float)pbc_.rHigh_.y_, (float)pbc_.rLow_.z_);
        glEnd();
    }

    for(double y=pbc_.rLow_.y_; y<=(pbc_.rHigh_.y_+stepY/2.0); y+=stepY)
    {
        glBegin(GL_LINES);
        glVertex3f((float)pbc_.rLow_.x_, (float)y, (float)pbc_.rLow_.z_);
        glVertex3f((float)pbc_.rLow_.x_, (float)y, (float)pbc_.rHigh_.z_);
        glEnd();
        
        glBegin(GL_LINES);
        glVertex3f((float)pbc_.rLow_.x_, (float)y, (float)pbc_.rLow_.z_);
        glVertex3f((float)pbc_.rHigh_.x_, (float)y, (float)pbc_.rLow_.z_);
        glEnd();        
    }

    for(double z=pbc_.rLow_.z_; z<=(pbc_.rHigh_.z_ + stepZ/2.0); z+=stepZ)
    {
        glBegin(GL_LINES);
        glVertex3f((float)pbc_.rLow_.x_, (float)pbc_.rLow_.y_, (float)z);
        glVertex3f((float)pbc_.rLow_.x_, (float)pbc_.rHigh_.y_, (float)z);
        glEnd();
        
        glBegin(GL_LINES);
        glVertex3f((float)pbc_.rLow_.x_, (float)pbc_.rLow_.y_, (float)z);
        glVertex3f((float)pbc_.rHigh_.x_, (float)pbc_.rLow_.y_, (float)z);
        glEnd();
    }


    textGap = maxSize*0.3;

    // Draw X-axis
    glPushMatrix();
    glTranslatef(float(pbc_.rHigh_.x_ + textGap), (float)pbc_.rLow_.y_, (float)pbc_.rLow_.z_);
    glRotatef(-45.0f, 1.0f, 0.0f, 0.0f);
    glScalef(float(maxSize*0.003), float(maxSize*0.003), float(maxSize*0.003));
    
    glColor3f(0.6f, 0.4f, 0.4f);
    glutStrokeCharacter(GLUT_STROKE_ROMAN, '-');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, '>');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, ' ');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, 'x');
    glPopMatrix();

    // Draw Y-axis
    glPushMatrix();
    glTranslatef((float)pbc_.rLow_.x_, float(pbc_.rHigh_.y_ + textGap), (float)pbc_.rLow_.z_);
    glRotatef(90.0f, 0.0f, 0.0f, 1.0f);
    glRotatef(45.0f, 1.0f, 0.0f, 0.0f);
    glScalef(float(maxSize*0.003), float(maxSize*0.003), float(maxSize*0.003));
    
    glColor3f(0.4f, 0.6f, 0.4f);
    glutStrokeCharacter(GLUT_STROKE_ROMAN, '-');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, '>');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, ' ');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, 'y');
    glPopMatrix();    

    // Draw Z-axis
    glPushMatrix();
    glTranslatef((float)pbc_.rLow_.x_, (float)pbc_.rLow_.y_, float(pbc_.rHigh_.z_ + textGap));
    glRotatef(-90.0f, 0.0f, 1.0f, 0.0f);
    glRotatef(225.0f, 1.0f, 0.0f, 0.0f);
    glScalef(float(maxSize*0.003), float(maxSize*0.003), float(maxSize*0.003));
    
    glColor3f(0.4f, 0.4f, 0.6f);
    glutStrokeCharacter(GLUT_STROKE_ROMAN, '-');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, '>');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, ' ');
    glutStrokeCharacter(GLUT_STROKE_ROMAN, 'z');
    glPopMatrix();    
    
    glEnable(GL_LIGHTING);
}

void C3DView::drawMolecules(ESelModes mode)
{
    C3DVector bond, end, emptySelColor;
    CAtom* A1;
    CAtom* A2;
    double r=0.0, g=0.0, b=0.0;
    double bondRadius;
    int selCount = 0;
    
    if(atoms_ && currentFrame_ && defaultAtProp_)
    {
        for(int i=0; i<atoms_->size(); i++)
        {
            std::string ID = (*atoms_)[i]->getID();
            defaultAtProp_->getCPKColor(ID.data(), r, g, b);

            if(mode == selmodeAtoms) glLoadName(i);
            drawAtom(currFrmAtVec(*(*atoms_)[i]), 0.3, (float)r, (float)g, (float)b);
            
            if(mode != selmodeAtoms)
            {
                A1 = (*atoms_)[i].get();
                for(int j=0; j<A1->getNumBonds(); j++)
                {
                    A2 = A1->getBondDest(j);
                 
                    bool bondAcrossPBC = false;
                    bond = currFrmAtVec(*A2) - currFrmAtVec(*A1);
                    if((fabs(bond.x_) > (pbc_.getWidthX()/2.0)) ||
                       (fabs(bond.y_) > (pbc_.getWidthY()/2.0)) ||
                       (fabs(bond.z_) > (pbc_.getWidthZ()/2.0)))    { bondRadius = 0.05; bondAcrossPBC = true; }
                    else                                              bondRadius = 0.2;
                    bond*= 0.5;
                    
                    end = currFrmAtVec(*A1) + bond;
                    
                    if(viewBondsAcrossPBC_ || !bondAcrossPBC)
                        drawBond(currFrmAtVec(*A1), end, bondRadius, (float)r, (float)g, (float)b);
                }
            }
        }

        if(mode != selmodeAtoms)
        {
            for(int i=0; i<atoms_->size(); i++)
            {            
                if((*atoms_)[i]->isSelected())
                {
                    if(selCount >= numSelRotColors_) selCount = 0;
                    drawAtom(currFrmAtVec(*(*atoms_)[i]), 0.3, (float)r, (float)g, (float)b, true, selColorRot_[selCount]);
                    selCount++;
                }
            }
        }
    }
}

void C3DView::drawGLObjects(ESelModes mode)
{
    if(mode == selmodeAtoms) return;
    
    if(glObjects_)
    {
        for(int i=0; i<glObjects_->size(); i++)
        {
            if((*glObjects_)[i]->getType() == CGLObject::objLine)
            {
                CGLObjectLine* obj = (CGLObjectLine*)(*glObjects_)[i].get();
                
                if(!obj) continue;
                
                glDisable(GL_LIGHTING);
                glColor4f(obj->color_[0], obj->color_[1], obj->color_[2], 0.0);

                glBegin(GL_LINES);
                glVertex3f((float)obj->p1_.x_, (float)obj->p1_.y_, (float)obj->p1_.z_);
                glVertex3f((float)obj->p2_.x_, (float)obj->p2_.y_, (float)obj->p2_.z_);
                glEnd();

                glEnable(GL_LIGHTING);
            }
        }
    }
}

void C3DView::drawSelectionHUD(ESelModes mode)
{
    C3DVector r;
    double hudXPos=-0.9, hydYPos=0.9;
    double X=hudXPos, Y=hydYPos;
    double lastXRasterPos, largestXRasterPos = -1.0;
    bool hasSel=false;
    std::string text;
    int selCount = 0;
    
    if(mode == selmodeAtoms) return;

    // Draw HUD text
    if(atoms_ && currentFrame_)
    {
        for(int i=0; i<atoms_->size(); i++)
        {
            if((*atoms_)[i]->isSelected() && (selCount < numSelRotColors_))
            {
                hasSel = true;

                std::string ID = (*atoms_)[i]->getID();
                text.resize(ID.length() + 512);
                r = currFrmAtVec(*(*atoms_)[i]);
                sprintf((char*)text.data(), "Atom sel [AtInd=%i, MolInd=%i, ID=%s (%s), pos=(%g,%g,%g), q=%g, m=%g]",
                        i, (*atoms_)[i]->getMolIndex(), ID.data(), defaultAtProp_ ? defaultAtProp_->getRecognizedAtomType(ID.data()).data() : "??",
                        r.x_, r.y_, r.z_, (*atoms_)[i]->Q_, (*atoms_)[i]->m_);
                lastXRasterPos = drawBitmapText(text.data(), GLUT_BITMAP_HELVETICA_12, X, Y, selColorRot_[selCount]);

                if(lastXRasterPos > largestXRasterPos) largestXRasterPos = lastXRasterPos;
                Y-=0.04;
                
                selCount++;
            }
        }

        // Draw current frame
        if((atoms_->size() > 0) && ((*atoms_)[0]->r_.size() > 1))
        {
            text.resize(64);
            sprintf((char*)text.data(), "Frame %i / %i", *currentFrame_, (int)(*atoms_)[0]->r_.size());
            drawBitmapText(text.data(), GLUT_BITMAP_HELVETICA_12, 0.65, -0.9, C3DVector(0.7, 0.7, 0.7));
        }
    }
    
    // Draw top and bottom line of HUD in
    // same projection and view matrix as
    // text
    if(hasSel)
    {
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glDisable(GL_DEPTH_TEST); 
        glDisable(GL_LIGHTING);
        glColor4f(0.4f, 0.4f, 0.4f, 1.0f);
        
        glBegin(GL_LINES);
        glVertex3d(hudXPos-0.03, hydYPos+0.04, 0.0);
        glVertex3d(largestXRasterPos+0.03, hydYPos+0.04, 0.0);
        glVertex3d(largestXRasterPos+0.03, Y+0.02, 0.0);
        glVertex3d(hudXPos-0.03, Y+0.02, 0.0);
        glEnd();
        
        glEnable(GL_LIGHTING);
        glEnable(GL_DEPTH_TEST);
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
    }
}

double C3DView::calcVDWDensity(double x, double y, double z)
{
    double alpha, dR2, density;
    C3DVector r(x, y, z);
    const double cutoffR2 = 9.0;
    const double scale = 0.5;

    density = 0.0;
    for(int i=0; i<atoms_->size(); i++)
    {
        if(!(*atoms_)[i]) continue;
        C3DVector R = r - currFrmAtVec(*(*atoms_)[i]);
        dR2 = R.norm2();
        if(dR2 > cutoffR2) continue;
        
        alpha = (*atoms_)[i]->sigma_ * scale;
        density+= expLookup_.exp(-dR2 / (2.0 * alpha));
    }
    
    return density;
}

void C3DView::drawIsoSurfaces(ESelModes mode)
{
    if(!viewIsoSurface_) return;
    
    C3DRect pbc(C3DVector((double)FLT_MAX, (double)FLT_MAX, (double)FLT_MAX), C3DVector(-(double)FLT_MAX, -(double)FLT_MAX, -(double)FLT_MAX));
    const double isoLevel = 0.606531; // e^{-1/2}
      
    if(mode == selmodeAtoms) return;
    if(!atoms_ || !currentFrame_) return;
    
    // Detect PBC of atom selection
    C3DVector r;
    for(int i=0; i<atoms_->size(); i++)
    {
        if(!(*atoms_)[i]) continue;
        
        r = currFrmAtVec(*(*atoms_)[i]);
        
        if(r.x_ > pbc.rHigh_.x_) pbc.rHigh_.x_ = r.x_;
        if(r.y_ > pbc.rHigh_.y_) pbc.rHigh_.y_ = r.y_;
        if(r.z_ > pbc.rHigh_.z_) pbc.rHigh_.z_ = r.z_;
        
        if(r.x_ < pbc.rLow_.x_) pbc.rLow_.x_ = r.x_;
        if(r.y_ < pbc.rLow_.y_) pbc.rLow_.y_ = r.y_;
        if(r.z_ < pbc.rLow_.z_) pbc.rLow_.z_ = r.z_;
    }
    
    // Add additional PBC space to take into account a finite
    // van der Waals R of the atoms located at the boundary
    // (assume maximum R of 3Å)
    pbc.rLow_.x_-= 3.0; pbc.rHigh_.x_+= 3.0;
    pbc.rLow_.y_-= 3.0; pbc.rHigh_.y_+= 3.0;
    pbc.rLow_.z_-= 3.0; pbc.rHigh_.z_+= 3.0;
    
    // Determine the Voxel dimensions
    double voxelX = pbc.getWidthX() / 32.0;
    double voxelY = pbc.getWidthY()/ 32.0;
    double voxelZ = pbc.getWidthZ() / 32.0;
    CMarchingCubesTriangle triangles[16];
    
    // Scan Voxel volume to determine iso-surfaces.
    // A density is assigned using a gaussian approx.,
    // see [EuroVis 2012 Short Papers, 14, 67 (2012)].    
    CMarchingCubesVoxel voxel;
    for(double x=pbc.rLow_.x_; x<=pbc.rHigh_.x_; x+=voxelX)
    {
        for(double y=pbc.rLow_.y_; y<=pbc.rHigh_.y_; y+=voxelY)
        {
            for(double z=pbc.rLow_.z_; z<=pbc.rHigh_.z_; z+=voxelZ)
            {
                voxel.val_[0] = calcVDWDensity(x, y, z);
                voxel.R_[0] = C3DVector(x, y, z);

                voxel.val_[1] = calcVDWDensity(x+voxelX, y, z);
                voxel.R_[1] = C3DVector(x+voxelX, y, z);

                voxel.val_[2] = calcVDWDensity(x+voxelX, y+voxelY, z);
                voxel.R_[2] = C3DVector(x+voxelX, y+voxelY, z);

                voxel.val_[3] = calcVDWDensity(x, y+voxelY, z);
                voxel.R_[3] = C3DVector(x, y+voxelY, z);

                voxel.val_[4] = calcVDWDensity(x, y, z+voxelZ);
                voxel.R_[4] = C3DVector(x, y, z+voxelZ);

                voxel.val_[5] = calcVDWDensity(x+voxelX, y, z+voxelZ);
                voxel.R_[5] = C3DVector(x+voxelX, y, z+voxelZ);
                
                voxel.val_[6] = calcVDWDensity(x+voxelX, y+voxelY, z+voxelZ);
                voxel.R_[6] = C3DVector(x+voxelX, y+voxelY, z+voxelZ);
                
                voxel.val_[7] = calcVDWDensity(x, y+voxelY, z+voxelZ);
                voxel.R_[7] = C3DVector(x, y+voxelY, z+voxelZ);


                int numTriangles = CMarchingCubes::calcTriangles(voxel, isoLevel, triangles);
                double eCoulomb;
                float color[4] = { 0.0f, 0.0f, 0.0f, 0.6f };
                
                
                CMolTwisterStateTools stateTools(nullptr, nullptr);
                for(int i=0; i<numTriangles; i++)
                {
                    glEnable(GL_BLEND);
                    glDepthMask(GL_FALSE);
                    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
                    
                    glBegin(GL_TRIANGLES);

                    eCoulomb = stateTools.measureCoulombPotential(triangles[i].R_[0], atoms_, *currentFrame_);
                    if(eCoulomb > 0.0) { color[0]=1.0f; color[1]=1.0f-(float)eCoulomb/0.001f; color[2]=1.0f-(float)eCoulomb/0.001f; }
                    else               { color[0]=1.0f+(float)eCoulomb/0.001f; color[1]=1.0f+(float)eCoulomb/0.001f; color[2]=1.0f; }
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
                    glNormal3f((float)triangles[i].N_[0].x_, (float)triangles[i].N_[0].y_, (float)triangles[i].N_[0].z_);
                    glVertex3f((float)triangles[i].R_[0].x_, (float)triangles[i].R_[0].y_, (float)triangles[i].R_[0].z_);
                    
                    eCoulomb = stateTools.measureCoulombPotential(triangles[i].R_[1], atoms_, *currentFrame_);
                    if(eCoulomb > 0.0) { color[0]=1.0f; color[1]=1.0f-(float)eCoulomb/0.001f; color[2]=1.0f-(float)eCoulomb/0.001f; }
                    else               { color[0]=1.0f+(float)eCoulomb/0.001f; color[1]=1.0f+(float)eCoulomb/0.001f; color[2]=1.0f; }
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
                    glNormal3f((float)triangles[i].N_[1].x_, (float)triangles[i].N_[1].y_, (float)triangles[i].N_[1].z_);
                    glVertex3f((float)triangles[i].R_[1].x_, (float)triangles[i].R_[1].y_, (float)triangles[i].R_[1].z_);

                    eCoulomb = stateTools.measureCoulombPotential(triangles[i].R_[2], atoms_, *currentFrame_);
                    if(eCoulomb > 0.0) { color[0]=1.0f; color[1]=1.0f-(float)eCoulomb/0.001f; color[2]=1.0f-(float)eCoulomb/0.001f; }
                    else               { color[0]=1.0f+(float)eCoulomb/0.001f; color[1]=1.0f+(float)eCoulomb/0.001f; color[2]=1.0f; }
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
                    glNormal3f((float)triangles[i].N_[2].x_, (float)triangles[i].N_[2].y_, (float)triangles[i].N_[2].z_);
                    glVertex3f((float)triangles[i].R_[2].x_, (float)triangles[i].R_[2].y_, (float)triangles[i].R_[2].z_);
                    
                    glEnd();

                    glDepthMask(GL_TRUE);
                    glDisable(GL_BLEND);
                }
            }
        }
    }
}

void C3DView::drawScene(ESelModes mode)
{
    C3DVector u, newPos;
      
    // Clear screen and load identity into transformation matrix
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // Rescale camera parameters based on zoom level and update camera
    u = camera_.lookAt_ - camera_.pos_;
    u*= camera_.zoomFactor_;
    
    newPos = camera_.lookAt_ - u;
    
    gluLookAt(newPos.x_, newPos.y_, newPos.z_,
              camera_.lookAt_.x_, camera_.lookAt_.y_, camera_.lookAt_.z_, 
              camera_.up_.x_, camera_.up_.y_, camera_.up_.z_);

    if(fogEnabled_)
    {
        glEnable(GL_FOG);
        C3DVector vLength = camera_.lookAt_ - newPos;
        float length = (float)vLength.norm();
        glFogf(GL_FOG_START, length / 2.0f);
        glFogf(GL_FOG_END, length);
    }
    else
    {
        glDisable(GL_FOG);
    }

    // Draw scene objects
    drawPBCGrid(mode);

    bool objectMoveMode = leftMButtonPressed_ || middleMButtonPressed_ || rightMButtonPressed_;
    bool movingAndLarge = objectMoveMode && atoms_ && atoms_->size() > numAtomsBeforeNoDraw_;
    if(!movingAndLarge)
    {
        drawMolecules(mode);
        drawGLObjects(mode);
        drawSelectionHUD(mode);
        drawIsoSurfaces(mode);
    }
}

void C3DView::show(std::vector<std::shared_ptr<CAtom>>* atoms, std::vector<std::shared_ptr<CGLObject>>* glObjects, int* currentFrame, CDefaultAtomicProperties* defAtProp)
{
    atoms_ = atoms;
    glObjects_ = glObjects;
    currentFrame_ = currentFrame;
    defaultAtProp_ = defAtProp;

    initOpenGL();
    initScene();
    update(true);
    
    glutMainLoop();
}

void C3DView::update(bool updateCameraPos)
{    
    double low, high;
    
    pbc_ = calcPBC();

    // Recalculate camera positions and frustum
    recalcFrustum(low, high);
    
    if(updateCameraPos)
    {
        camera_.pos_.x_ = pbc_.rHigh_.x_ + (high - low);
        camera_.pos_.y_ = (pbc_.rHigh_.y_ + pbc_.rLow_.y_) / 2.0;
        camera_.pos_.z_ = (pbc_.rHigh_.z_ + pbc_.rLow_.z_) / 2.0;
        
        camera_.lookAt_ = C3DVector((pbc_.rHigh_.x_ + pbc_.rLow_.x_) / 2.0, 
                                    (pbc_.rHigh_.y_ + pbc_.rLow_.y_) / 2.0,
                                    (pbc_.rHigh_.z_ + pbc_.rLow_.z_) / 2.0);
        camera_.up_ = C3DVector(0.0, 0.0, 1.0);
    }

    glutPostRedisplay();
}

void C3DView::recalcFrustum(double& frustLow, double& frustHigh, ESelModes mode)
{
    double low = pbc_.rLow_.x_;
    double high = pbc_.rHigh_.x_;
    double nearClip, farClip;
    double halfDist;
    double aspectRatio1, aspectRatio2;
    int w = lastWindowSize_.x_;
    int h = lastWindowSize_.y_;
      
    if(pbc_.rLow_.y_ < low) low = pbc_.rLow_.y_;
    if(pbc_.rLow_.z_ < low) low = pbc_.rLow_.z_;
    if(pbc_.rHigh_.y_ > high) high = pbc_.rHigh_.y_;
    if(pbc_.rHigh_.z_ > high) high = pbc_.rHigh_.z_;
    low*= 1.1;
    high*= 1.1;
    
    nearClip = 0.5*(high - low);
    farClip = 3.0*(high - low);
    
    nearClip*= camera_.zoomFactor_;
    farClip*= camera_.zoomFactor_;
    
    halfDist = 0.5*(high - low)*camera_.zoomFactor_;

    if(mode == selmodeNone)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
    }
    if(w <= h)
    {
        aspectRatio1 = 1.0;
        aspectRatio2 = ((w == 0) ? 1.0 : double(h) / double(w));
    }
    else
    {
        aspectRatio1 = ((w == 0) ? 1.0 : double(w) / double(h));
        aspectRatio2 = 1.0;
    }
    if(orthoView_) glOrtho(-halfDist*aspectRatio1, halfDist*aspectRatio1, -halfDist*aspectRatio2, halfDist*aspectRatio2, nearClip, farClip);
    else           glFrustum(-halfDist*aspectRatio1, halfDist*aspectRatio1, -halfDist*aspectRatio2, halfDist*aspectRatio2, nearClip, farClip);
    
    if(mode == selmodeNone)
    {
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
    
    frustLow = low;
    frustHigh = high;
}

void C3DView::onRender()
{
    drawScene(selmodeNone);
    
    glutSwapBuffers();
}

void C3DView::onReshape(int w, int h)
{
    double  low, high;
      
    lastWindowSize_ = CScreenCoord(w, h);
    glViewport(0, 0, w, h);
    
    recalcFrustum(low, high);
}

void C3DView::onTimer(int)
{
    if(updateRequested_ == 1)
    {
        update(false);
        updateRequested_ = 0;
    }
    
    if(updateRequested_ == 2)
    {
        update(true);
        updateRequested_ = 0;
    }
    
    if(fullscreenRequested_ == 1)
    {
        glutFullScreen();
        fullscreenRequested_ = 0;
    }
    
    if(fullscreenRequested_ == 2)
    {
        int W = glutGet(GLUT_SCREEN_WIDTH);
        int H = glutGet(GLUT_SCREEN_HEIGHT);
        glutPositionWindow(W/2, 0);
        glutReshapeWindow(W/2, H);
        fullscreenRequested_ = 0;
    }
    
    if(requestQuit_)
    {
        glutDestroyWindow(glutGetWindow());
        exit(0);
    }

    glutTimerFunc(100, onTimer, 0);
}

void C3DView::onMouseClick(int button, int state, int x, int y)
{
    int modifier = glutGetModifiers();

    if(state == GLUT_DOWN) coordLastClick_ = CScreenCoord(x, y);

    if(modifier == GLUT_ACTIVE_SHIFT)
    {
        if(button == GLUT_LEFT_BUTTON)
        {
            if(state == GLUT_UP)
            {
                pickAtoms(x, y);
            }
        }
        if(button == GLUT_RIGHT_BUTTON)
        {
            if(state == GLUT_UP)
            {
                if(atoms_)
                {
                    for(int i=0; i<atoms_->size(); i++)
                    {
                        (*atoms_)[i]->select(false);
                    }
                }

                glutPostRedisplay();
            }
        }
    }
    
    else
    {
        if(button == GLUT_LEFT_BUTTON)
        {
            if(state == GLUT_DOWN) leftMButtonPressed_ = true;
            else                   leftMButtonPressed_ = false;

            if(!leftMButtonPressed_) glutPostRedisplay();
        }
        
        if(button == GLUT_MIDDLE_BUTTON)
        {
            if(state == GLUT_DOWN) middleMButtonPressed_ = true;
            else                   middleMButtonPressed_ = false;

            if(!middleMButtonPressed_) glutPostRedisplay();
        }
        
        if(button == GLUT_RIGHT_BUTTON)
        {
            if(state == GLUT_DOWN) rightMButtonPressed_ = true;
            else                   rightMButtonPressed_ = false;

            if(!rightMButtonPressed_) glutPostRedisplay();
        }
    }
}

void C3DView::onMouseMove(int x, int y)
{
    if(leftMButtonPressed_)
    {
        C3DVector L = camera_.pos_ - camera_.lookAt_;
        double absL = L.norm();
        double Dx = coordLastClick_.x_ - x;
        double Dy = y - coordLastClick_.y_;
        double low, high;
        double theta, phi;

        if((lastWindowSize_.x_ != 0.0) && (absL != 0.0))
        {
            if(fabs(Dx) > fabs(Dy))
            {
                theta = 2.0 * M_PI * (Dx / lastWindowSize_.x_);
                camera_.pos_ = camera_.lookAt_ + L*cos(theta) + camera_.up_.cross(L)*sin(theta);
            }
            else
            {
                phi = 2.0 * M_PI * (Dy / lastWindowSize_.y_);
                camera_.pos_ = camera_.lookAt_ + L*cos(phi) + camera_.up_*(absL*sin(phi));
                camera_.up_ = camera_.up_*cos(phi) - L*(sin(phi) / absL);
            }

            recalcFrustum(low, high);
            glutPostRedisplay();
        }

        coordLastClick_ = CScreenCoord(x, y);
    }
    
    else if(middleMButtonPressed_)
    {
        double D = (coordLastClick_.y_ - y) + (x - coordLastClick_.x_);
        double low, high;
        
        if(lastWindowSize_.y_ != 0.0)
        {
            camera_.zoomFactor_-= D / lastWindowSize_.y_;
            if(camera_.zoomFactor_ > 10.0) camera_.zoomFactor_ = 10.0;
            if(camera_.zoomFactor_ < 0.1) camera_.zoomFactor_ = 0.1;
            recalcFrustum(low, high);
            glutPostRedisplay();
        }
        
        coordLastClick_ = CScreenCoord(x, y);
    }

    else if(rightMButtonPressed_)
    {
        C3DVector L, u, M;
        double Dx = coordLastClick_.x_ - x;
        double Dy = y - coordLastClick_.y_;
        double low, high, absL;

        if((lastWindowSize_.x_ != 0.0) && (lastWindowSize_.y_ != 0.0))
        {
            L = camera_.lookAt_ - camera_.pos_;
            absL = L.norm();
            
            if(absL != 0.0)
            {
                u = L.cross(camera_.up_)*(1.0 / absL);
                
                Dx = (Dx / lastWindowSize_.x_)*pbc_.getLargestWidth();
                Dy = (Dy / lastWindowSize_.y_)*pbc_.getLargestWidth();
                
                M = u*Dx + camera_.up_*Dy;
                camera_.pos_ = camera_.pos_ + M;
                camera_.lookAt_ = camera_.lookAt_ + M;

                recalcFrustum(low, high);
                glutPostRedisplay();
            }
        }

        coordLastClick_ = CScreenCoord(x, y);
    }
}

void C3DView::onKeyboard(unsigned char key, int, int)
{
    int mod = glutGetModifiers();
    
    if(key == 27) // Esc-key
    {
        fullscreenRequested_ = 2;
        glutPostRedisplay();
    }
    
    if(key == 'f')
    {
        if(mod & GLUT_ACTIVE_ALT)
        {
            fullscreenRequested_ = 1;
            glutPostRedisplay();
        }
    }

    if(key == 'o')
    {
        if(mod & GLUT_ACTIVE_ALT)
        {
            if(orthoView_) orthoView_ = false;
            else           orthoView_ = true;

            updateRequested_ = 1;
            glutPostRedisplay();
        }
    }

    if(key == 'a')
    {
        if(mod & GLUT_ACTIVE_ALT)
        {
            if(viewAxes_) viewAxes_ = false;
            else          viewAxes_ = true;
            
            glutPostRedisplay();
        }
    }

    if(key == 'i')
    {
        if(mod & GLUT_ACTIVE_ALT)
        {
            if(viewIsoSurface_) viewIsoSurface_ = false;
            else                viewIsoSurface_ = true;
            
            glutPostRedisplay();
        }
    }
    
    if(key == 'p')
    {
        if(mod & GLUT_ACTIVE_ALT)
        {
            if(viewBondsAcrossPBC_) viewBondsAcrossPBC_ = false;
            else                    viewBondsAcrossPBC_ = true;
            
            glutPostRedisplay();
        }
    }
}

void C3DView::onSpecialFunc(int key, int, int)
{
    if(key == GLUT_KEY_RIGHT)
    {
        if(currentFrame_ && atoms_)
        {
            (*currentFrame_)++;
            if((atoms_->size() > 0) && ((*currentFrame_) >= (*atoms_)[0]->r_.size()))
                (*currentFrame_) = 0;

            update(false);
            glutPostRedisplay();
        }
    }

    if(key == GLUT_KEY_LEFT)
    {
        if(currentFrame_ && atoms_)
        {
            (*currentFrame_)--;
            if((atoms_->size() > 0) && ((*currentFrame_) < 0))
                (*currentFrame_) = (int)(*atoms_)[0]->r_.size()-1;
            
            update(false);
            glutPostRedisplay();
        }
    }
}

void C3DView::pickAtoms(int x, int y)
{
    C3DVector u, newPos;
    std::vector<GLuint> selectBuffer(atoms_->size());
    GLuint* ptr = selectBuffer.data();
    double frustLow, frustHigh;
    GLint hits;
    GLint viewport[4];
    

    // Enter selection mode
    glSelectBuffer(((int)atoms_->size() + 1)*sizeof(GLuint), selectBuffer.data());
    glRenderMode(GL_SELECT);
    
    // Switch to projection matrix mode and reset matrix
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    // Retrieve viewport dimensions (x, y, width, height) and
    // use this info to select a projection that only covers a
    // small 5x5 area of the view.
    glGetIntegerv(GL_VIEWPORT, viewport);
    gluPickMatrix(x, viewport[3]-y, 5.0, 5.0, viewport);
    recalcFrustum(frustLow, frustHigh, selmodeAtoms);
    

    // Draw objects that can be selected in selection mode
    glMatrixMode(GL_MODELVIEW);
    glInitNames();
    glPushName(0);
    drawScene(selmodeAtoms);
    
    
    // Restore old projection matrix and flush modelview pipeline
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glFlush();
    

    
    // Switch back to render mode (from sel. mode)
    // and search for object in selection that is
    // closest to the viewport (i.e. smalledt depth
    // value in the z-buffer)
    float fZ1, fZ2, minZ=0.0f;
    int selAtIndex=-1;
    
    hits = glRenderMode(GL_RENDER);
    for(int i=0; i<hits; i++)
    {
        int numNames, currSel=0;
        
        numNames = *ptr;
        
        ptr++;
        fZ1 = float(*ptr) / float(0xFFFFFFFF);
        ptr++;
        fZ2 = float(*ptr) / float(0xFFFFFFFF);
        
        ptr++;
        for(int j=0; j<numNames; j++)
        {
            currSel = *ptr;
            ptr++;
        }
        
        if(i == 0) { minZ = fZ1; selAtIndex = currSel; }
        if(fZ1 < minZ) { minZ = fZ1; selAtIndex = currSel; }
    }
    
    if(atoms_)
    {
        if((selAtIndex != -1) && (selAtIndex < atoms_->size()))
        {
            if((*atoms_)[selAtIndex]->isSelected())
                (*atoms_)[selAtIndex]->select(false);
            else
                (*atoms_)[selAtIndex]->select(true);
        }
    }
    
    // Redisplay with changes
    glutPostRedisplay();
}
