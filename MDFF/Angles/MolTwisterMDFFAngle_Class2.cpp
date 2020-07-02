#include "MolTwisterMDFFAngle_Class2.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

void CMDFFAngle_Class2::serialize(std::stringstream& io, bool saveToStream)
{
    CMDFFAngle::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << Ea_.theta0_;
        io << Ea_.K2_;
        io << Ea_.K3_;
        io << Ea_.K4_;

        io << Ebb_.M_;
        io << Ebb_.r1_;
        io << Ebb_.r2_;

        io << Eba_.N1_;
        io << Eba_.N2_;
        io << Eba_.r1_;
        io << Eba_.r2_;

        io << toRad_;
    }
    else
    {
        io >> Ea_.theta0_;
        io >> Ea_.K2_;
        io >> Ea_.K3_;
        io >> Ea_.K4_;

        io >> Ebb_.M_;
        io >> Ebb_.r1_;
        io >> Ebb_.r2_;

        io >> Eba_.N1_;
        io >> Eba_.N2_;
        io >> Eba_.r1_;
        io >> Eba_.r2_;

        io >> toRad_;
    }
}

size_t CMDFFAngle_Class2::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;

    // Read the Ea coefficients
    Ea_.theta0_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    Ea_.K2_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    Ea_.K3_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    Ea_.K4_ = atof(CASCIIUtility::getArg(arguments, arg++).data());

    // Read the Ebb coefficients
    Ebb_.M_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    Ebb_.r1_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    Ebb_.r2_ = atof(CASCIIUtility::getArg(arguments, arg++).data());

    // Read the Eba coefficients
    Eba_.N1_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    Eba_.N2_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    Eba_.r1_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    Eba_.r2_ = atof(CASCIIUtility::getArg(arguments, arg++).data());

    return arg;
}

std::string CMDFFAngle_Class2::getArguments(bool convertToCal) const
{
    char string[512];
    
    sprintf(string, "Ea:{theta0:%.4f K2:%.4f K3:%.4f K4:%.4f} Ebb:{M:%.4f r1:%.4f r2:%.4f} Eba:{N1:%.4f N2:%.4f r1:%.4f r2:%.4f}",
            Ea_.theta0_, J2cal(Ea_.K2_, convertToCal), J2cal(Ea_.K3_, convertToCal), J2cal(Ea_.K4_, convertToCal),
            J2cal(Ebb_.M_, convertToCal), Ebb_.r1_, Ebb_.r2_,
            J2cal(Eba_.N1_, convertToCal), J2cal(Eba_.N2_, convertToCal), Eba_.r1_, Eba_.r2_);
    
    return string;
}

std::string CMDFFAngle_Class2::getArgumentsLaTeX(bool convertToCal) const
{
    char string[512];
    
    sprintf(string, "$E_a$:\\{$\\theta_0$=%g, $K_2$=%g, $K_3$=%g, $K_4$=%g\\}, $E_{bb}$:\\{$M$=%g, $r_1$=%g, $r_2$=%g\\}, $E_{ba}$:\\{$N_1$=%g, $N_2$=%g, $r_1$=%g, $r_2$=%g\\}",
            Ea_.theta0_, J2cal(Ea_.K2_, convertToCal), J2cal(Ea_.K3_, convertToCal), J2cal(Ea_.K4_, convertToCal),
            J2cal(Ebb_.M_, convertToCal), Ebb_.r1_, Ebb_.r2_,
            J2cal(Eba_.N1_, convertToCal), J2cal(Eba_.N2_, convertToCal), Eba_.r1_, Eba_.r2_);
    
    return string;
}

std::string CMDFFAngle_Class2::getLammpsDef(int index, bool convertToCal) const
{
    char line[1024];
    
    if(index == 0)
    {
        sprintf(line, "class2 %.6f %.6f %.6f %.6f     # %s %s %s", Ea_.theta0_, J2cal(Ea_.K2_, convertToCal), J2cal(Ea_.K3_, convertToCal), J2cal(Ea_.K4_, convertToCal),
                atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), atomNamesToBond_[2].data());
    }
    if(index == 1)
    {
        sprintf(line, "class2 bb %.6f %.6f %.6f     # %s %s %s", J2cal(Ebb_.M_, convertToCal), Ebb_.r1_, Ebb_.r2_,
                atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), atomNamesToBond_[2].data());
    }
    if(index == 2)
    {
        sprintf(line, "class2 ba %.6f %.6f %.6f %.6f     # %s %s %s", J2cal(Eba_.N1_, convertToCal), J2cal(Eba_.N2_, convertToCal), Eba_.r1_, Eba_.r2_,
                atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), atomNamesToBond_[2].data());
    }

    return line;
}

std::string CMDFFAngle_Class2::getHoomdBlueDef(int, int) const
{
    return "Not available!";
}

double CMDFFAngle_Class2::calcPotentialHarm(C3DVector r1, C3DVector r2, C3DVector r3, double k, double theta0, int degree) const
{
    C3DVector   r21 = r1 - r2;
    C3DVector   r23 = r3 - r2;
    double      R21 = r21.norm();
    double      R23 = r23.norm();
    double      R21Inv = (R21 == 0.0) ? 1.0 / 1E-10 : 1.0 / R21;
    double      R23Inv = (R23 == 0.0) ? 1.0 / 1E-10 : 1.0 / R23;
    double      cosTheta = (r21*r23) * (R21Inv*R23Inv);
    double      Theta = acos(cosTheta);
    double      Theta0 = theta0 * toRad_;
    double      DeltaTheta = Theta - Theta0;
    
    double DeltaN = 1.0;
    for(int i=0; i<degree; i++)
        DeltaN*= DeltaTheta;
    
    return k * DeltaN;
}

void CMDFFAngle_Class2::calcForcesHarm(C3DVector r1, C3DVector r2, C3DVector r3, double k, double theta0, int degree, C3DVector& f1, C3DVector& f2, C3DVector& f3) const
{
    C3DVector   r21 = r1 - r2;
    C3DVector   r23 = r3 - r2;
    double      R21 = r21.norm();
    double      R23 = r23.norm();
    double      R21Inv = (R21 == 0.0) ? 1.0 / 1E-10 : 1.0 / R21;
    double      R23Inv = (R23 == 0.0) ? 1.0 / 1E-10 : 1.0 / R23;
    double      cosTheta = (r21*r23) * (R21Inv*R23Inv);
    double      Theta = acos(cosTheta);
    double      Theta0 = theta0 * toRad_;
    double      deltaTheta = Theta - Theta0;
    C3DVector   dTheta_dr3 = calcAngularForceCoeffs13(r1, r2, r3);
    C3DVector   dTheta_dr2 = calcAngularForceCoeffs2(r1, r2, r3);
    C3DVector   dTheta_dr1 = calcAngularForceCoeffs13(r3, r2, r1);
    double      dU_dTheta = -double(degree)*k;

    for(int i=0; i<(degree-1); i++)
        dU_dTheta*= deltaTheta;
    
    f3 = dTheta_dr3*dU_dTheta;
    f2 = dTheta_dr2*dU_dTheta;
    f1 = dTheta_dr1*dU_dTheta;
}

double CMDFFAngle_Class2::calcPotentialCartesian(C3DVector r1, C3DVector r2, C3DVector r3, double M, double ra, double rb) const
{
    C3DVector   r12 = r2 - r1;
    C3DVector   r23 = r3 - r2;
    double      R12 = r12.norm();
    double      R23 = r23.norm();
    
    return M * (R12 - ra) * (R23 - rb);
}

void CMDFFAngle_Class2::calcForcesCartesian(C3DVector r1, C3DVector r2, C3DVector r3, double M, double ra, double rb, C3DVector& f1, C3DVector& f2, C3DVector& f3) const
{
    C3DVector   r12 = r2 - r1;
    C3DVector   r23 = r3 - r2;
    double      R12 = r12.norm();
    double      R23 = r23.norm();
    double      RInv12 = (R12 == 0.0) ? 1.0 / 1E-10 : 1.0 / R12;
    double      RInv23 = (R23 == 0.0) ? 1.0 / 1E-10 : 1.0 / R23;
    double      R12_ra = R12 - ra;
    double      R23_rb = R23 - rb;
    C3DVector   r12xRInv12 = r12*RInv12;
    C3DVector   r23xRInv23 = r23*RInv23;
    
    f1 = r12xRInv12*(M*R23_rb);
    f3 = r23xRInv23*((-1.0)*M*R12_ra);
    f2 = (f1 + f3)*(-1.0);
}

double CMDFFAngle_Class2::calcPotentialMix(C3DVector r1, C3DVector r2, C3DVector r3, double N1, double N2, double ra, double rb, double theta0) const
{
    C3DVector   r21 = r1 - r2;
    C3DVector   r23 = r3 - r2;
    double      R21 = r21.norm();
    double      R23 = r23.norm();
    double      R21Inv = (R21 == 0.0) ? 1.0 / 1E-10 : 1.0 / R21;
    double      R23Inv = (R23 == 0.0) ? 1.0 / 1E-10 : 1.0 / R23;
    double      cosTheta = (r21*r23) * (R21Inv*R23Inv);
    double      Theta = acos(cosTheta);
    double      Theta0 = theta0 * toRad_;
    double      deltaTheta = Theta - Theta0;
    
    return (N1*(R21 - ra) + N2*(R23 - rb)) * deltaTheta;
}

void CMDFFAngle_Class2::calcForcesMix(C3DVector r1, C3DVector r2, C3DVector r3, double N1, double N2, double ra, double rb, double theta0, C3DVector& f1, C3DVector& f2, C3DVector& f3) const
{
    C3DVector   r12 = r2 - r1;
    C3DVector   r23 = r3 - r2;
    C3DVector   r21 = r12*(-1.0);
    double      R12 = r12.norm();
    double      R23 = r23.norm();
    double      RInv12 = (R12 == 0.0) ? 1.0 / 1E-10 : 1.0 / R12;
    double      RInv23 = (R23 == 0.0) ? 1.0 / 1E-10 : 1.0 / R23;
    double      coeff = N1*(R12 - ra) + N2*(R23 - rb);
    double      cosTheta = (r21*r23) * (RInv12*RInv23);
    double      Theta = acos(cosTheta);
    double      Theta0 = theta0 * toRad_;
    double      deltaTheta = Theta - Theta0;
    C3DVector   r12xRInv12 = r12*RInv12;
    C3DVector   r23xRInv23 = r23*RInv23;
    C3DVector   dTheta_dr1 = calcAngularForceCoeffs13(r3, r2, r1);
    C3DVector   dTheta_dr2 = calcAngularForceCoeffs2(r1, r2, r3);
    C3DVector   dTheta_dr3 = calcAngularForceCoeffs13(r1, r2, r3);
    
    f1 = r12xRInv12*(N1*deltaTheta) - dTheta_dr1*coeff;
    f2 = r12xRInv12*((-1.0)*N1*deltaTheta) + r23xRInv23*(N2*deltaTheta) - dTheta_dr2*coeff;
    f3 = r23xRInv23*((-1.0)*N2*deltaTheta) - dTheta_dr3*coeff;
}

double CMDFFAngle_Class2::calcPotential(C3DVector r1, C3DVector r2, C3DVector r3) const
{
    double U = 0.0;
    
    U+= calcPotentialHarm(r1, r2, r3, Ea_.K2_, Ea_.theta0_, 2);
    U+= calcPotentialHarm(r1, r2, r3, Ea_.K3_, Ea_.theta0_, 3);
    U+= calcPotentialHarm(r1, r2, r3, Ea_.K4_, Ea_.theta0_, 4);
    
    U+= calcPotentialCartesian(r1, r2, r3, Ebb_.M_, Ebb_.r1_, Ebb_.r2_);

    U+= calcPotentialMix(r1, r2, r3, Eba_.N1_, Eba_.N2_, Eba_.r1_, Eba_.r2_, Ea_.theta0_);
    
    return U;
}

void CMDFFAngle_Class2::calcForces(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector& f1, C3DVector& f2, C3DVector& f3) const
{
    C3DVector F1, F2, F3;

    f1.set(0.0, 0.0, 0.0);
    f2.set(0.0, 0.0, 0.0);
    f3.set(0.0, 0.0, 0.0);

    calcForcesHarm(r1, r2, r3, Ea_.K2_, Ea_.theta0_, 2, F1, F2, F3);
    f1+= F1; f2+= F2; f3+= F3;
    calcForcesHarm(r1, r2, r3, Ea_.K3_, Ea_.theta0_, 3, F1, F2, F3);
    f1+= F1; f2+= F2; f3+= F3;
    calcForcesHarm(r1, r2, r3, Ea_.K4_, Ea_.theta0_, 4, F1, F2, F3);
    f1+= F1; f2+= F2; f3+= F3;

    calcForcesCartesian(r1, r2, r3, Ebb_.M_, Ebb_.r1_, Ebb_.r2_, F1, F2, F3);
    f1+= F1; f2+= F2; f3+= F3;

    calcForcesMix(r1, r2, r3, Eba_.N1_, Eba_.N2_, Eba_.r1_, Eba_.r2_, Ea_.theta0_, F1, F2, F3);
    f1+= F1; f2+= F2; f3+= F3;
}

END_CUDA_COMPATIBLE()
