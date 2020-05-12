#include "VelVerlet.h"
#include "Math.h"
#include <float.h>
/*
void CVelVerlet::Propagator(int N, int dim, double dt, double Lmax, vector<CParticle3D>& aParticles, vector<C3DVector> &aF, vector<C3DVector> &aFpi, bool bNPT)
{
    if(!bNPT)
    {
        double dt_2 = dt / 2.0;
        for(int k=0; k<N; k++)
        {
            double m_k = aParticles[k].m;
            aParticles[k].p+= aF[k]*dt_2;
            aParticles[k].x+= aParticles[k].p*(dt / m_k);
            aF[k] = CalcParticleForce(k, N, dim, Lmax, Lmax, Lmax, aParticles, aFpi[k]);
            aParticles[k].p+= aF[k]*dt_2;
        }
    }
    
    else
    {
        p_eps+= ((dt / 2.0) * G_eps(N, aParticles, aFpi));  // Step 3.1
        Prop_p(N, dt, aParticles, aF);                      // Step 3.2
        Prop_r(N, dt, aParticles, aF);                      // Step 3.3
        eps+= (dt * p_eps / W);                             // Step 3.4
        double Lm = pow(GetV(Lmax, true), 1.0/3.0);
        for(int k=0; k<N; k++)
            aF[k] = CalcParticleForce(k, N, dim, Lm, Lm, Lm, aParticles, aFpi[k]);
        Prop_p(N, dt, aParticles, aF);                      // Step 3.5
        p_eps+= ((dt / 2.0) * G_eps(N, aParticles, aFpi));  // Step 3.6
    }
}

void CVelVerlet::Prop_p(int N, double dt, vector<CParticle3D>& aParticles, vector<C3DVector> &aF)
{
    double u_eps = p_eps / W;
    double alpha = (1.0 + 1.0/double(N));
    double parm = alpha*u_eps*dt / 4.0;
    double coeff = CMt::Exp(-parm);
    double coeff2 = coeff*coeff;
    double coeff3 = CMt::SinhXoverX(parm); // sinh(parm) / parm;
    for(int k=0; k<N; k++)
    {
        aParticles[k].p = aParticles[k].p*coeff2 + aF[k]*((dt/2.0)*coeff*coeff3);
    }
}

void CVelVerlet::Prop_r(int N, double dt, vector<CParticle3D>& aParticles, vector<C3DVector> &aF)
{
    double u_eps = p_eps / W;
    double parm = u_eps*dt / 2.0;
    double coeff = CMt::Exp(parm);
    double coeff2 = coeff*coeff;
    double coeff3 = CMt::SinhXoverX(parm) * coeff;
    for(int k=0; k<N; k++)
    {
        C3DVector u_k = aParticles[k].p * (1.0 / aParticles[k].m);
        aParticles[k].x = aParticles[k].x*coeff2 + u_k*(dt*coeff3);
    }
}

C3DVector CVelVerlet::CalcParticleForce(int k, int N, int dim, double Lx, double Ly, double Lz, vector<CParticle3D>& aParticles, C3DVector& Fpi)
{
    C3DVector F;
    C3DVector PBCx = C3DVector( Lx, 0.0, 0.0);
    C3DVector PBCy = C3DVector(0.0,  Ly, 0.0);
    C3DVector PBCz = C3DVector(0.0, 0.0,  Lz);
    
    // Add external forces on particle (Fpi is force from primary image only)
    Fpi = C3DVector(0.0, 0.0, 0.0);
    Fpi+= aFExternal[k].CalcForce(aParticles[k].x);
    
    // Add non-bonded forces to particle, as well as
    // non-bonded forces from first PBC images
    for(int i=0; i<N; i++)
    {
        C3DVector r_k = aParticles[k].x;
        C3DVector r_i = aParticles[i].x;
        
        Fpi+= aFNonBonded[i][k].CalcForce(r_k, r_i);
        
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i + PBCx);
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i - PBCx);
        
        if(dim < 2) continue;
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i + PBCy);
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i - PBCy);
        
        if(dim < 3) continue;
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i + PBCz);
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i - PBCz);
    }
    F+= Fpi;
    
    // Add forces from harmonic bonds on particle k
    for(int j=0; j<aFHarmBond.size(); j++)
    {
        int iBondTo = -1;
        int iPotential = -1;
        if(aFHarmBond[j].m_i == k)
        {
            iBondTo = aFHarmBond[j].m_j;
            iPotential = j;
        }
        if(aFHarmBond[j].m_j == k)
        {
            iBondTo = aFHarmBond[j].m_i;
            iPotential = j;
        }
        if(iPotential != -1)
        {
            C3DVector r_k = aParticles[k].x;
            C3DVector r_i = aParticles[iBondTo].x;
            
            Fpi+= aFHarmBond[iPotential].CalcForce(r_k, r_i, Lx, Ly, Lz);
        }
    }
    F+= Fpi;
    
    if(fabs(F.x_) > Fcut) { F.x_ = ((F.x_ >= 0.0) ? 1.0 : -1.0) * Fcut; m_bCutF = true; }
    if(fabs(F.y_) > Fcut) { F.y_ = ((F.y_ >= 0.0) ? 1.0 : -1.0) * Fcut; m_bCutF = true; }
    if(fabs(F.z_) > Fcut) { F.z_ = ((F.z_ >= 0.0) ? 1.0 : -1.0) * Fcut; m_bCutF = true; }
    StoreMaxF(F);
    PrintDebugInfoAtCutForces(k, N, Lx, Ly, Lz, aParticles);
    
    return F;
}

double CVelVerlet::G_eps(int N, vector<CParticle3D>& aParticles, vector<C3DVector> &aF)
{
    double V = V0 * exp(3.0 * eps);
    
    double sum_p = 0.0;
    double sum_f = 0.0;
    for(int k=0; k<N; k++)
    {
        double p2 = aParticles[k].p * aParticles[k].p;
        sum_p+= (p2 / aParticles[k].m);
        sum_f+= (aF[k] * aParticles[k].x);
    }

    return (1.0 + 1.0 / double(N))*sum_p + sum_f - 3.0*P*V;
}

void CVelVerlet::SetRandMom(double tau)
{
    double          a = 2.0 / double(RAND_MAX);
    
    p_eps = (a*double(rand()) - 1.0) * (W / tau);
    if(p_eps == 0.0) p_eps = (W / tau);
}

double CVelVerlet::GetV(double Lmax, bool bNPT)
{
    if(!bNPT) return Lmax*Lmax*Lmax;
    
    return V0 * exp(3.0 * eps);
}

void CVelVerlet::StoreMaxF(C3DVector& F)
{
    double fAbs = F.Abs();
    
    if(fAbs > m_dLastFMax) m_dLastFMax = fAbs;
}

void CVelVerlet::PrintCutMsgAndReset()
{
    if(m_bCutF)
    {
//        COut::Printf("\t****************************************\r\n");
//        COut::Printf("\t* Warning! Forces were cut!!!          *\r\n");
//        COut::Printf("\t****************************************\r\n");
        m_dLastFMax = 0.0;
        m_bCutF = false;
    }
}

void CVelVerlet::PrintDebugInfoAtCutForces(int k, int N, double Lx, double Ly, double Lz, vector<CParticle3D>& aParticles)
{
    C3DVector F;
    C3DVector PBCx = C3DVector( Lx, 0.0, 0.0);
    C3DVector PBCy = C3DVector(0.0,  Ly, 0.0);
    C3DVector PBCz = C3DVector(0.0, 0.0,  Lz);

    if(m_bCutF)
    {
        C3DVector r_k = aParticles[k].x;
//        COut::Printf("\r\n\r\n------------------ Center -----------------------\r\n");
//        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x, r_k.y, r_k.z, F.x, F.y, F.z, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
//            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.Abs() > 10.0) ? "***" : "   ", i, r_i.x, r_i.y, r_i.z, F.x, F.y, F.z);
        }
        
//        COut::Printf("\r\n\r\n------------------ Plus PBCx -----------------------\r\n");
//        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x, r_k.y, r_k.z, F.x, F.y, F.z, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x + PBCx;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
//            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.Abs() > 10.0) ? "***" : "   ", i, r_i.x, r_i.y, r_i.z, F.x, F.y, F.z);
        }
//        COut::Printf("\r\n\r\n------------------ Minus PBCx -----------------------\r\n");
//        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x, r_k.y, r_k.z, F.x, F.y, F.z, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x - PBCx;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
//            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.Abs() > 10.0) ? "***" : "   ", i, r_i.x, r_i.y, r_i.z, F.x, F.y, F.z);
        }
//        COut::Printf("\r\n\r\n------------------ Plus PBCy -----------------------\r\n");
//        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x, r_k.y, r_k.z, F.x, F.y, F.z, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x + PBCy;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
//            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.Abs() > 10.0) ? "***" : "   ", i, r_i.x, r_i.y, r_i.z, F.x, F.y, F.z);
        }
//        COut::Printf("\r\n\r\n------------------ Minus PBCy -----------------------\r\n");
//        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x, r_k.y, r_k.z, F.x, F.y, F.z, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x - PBCy;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
//            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.Abs() > 10.0) ? "***" : "   ", i, r_i.x, r_i.y, r_i.z, F.x, F.y, F.z);
        }
//        COut::Printf("\r\n\r\n------------------ Plus PBCz -----------------------\r\n");
//        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x, r_k.y, r_k.z, F.x, F.y, F.z, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x + PBCz;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
//            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.Abs() > 10.0) ? "***" : "   ", i, r_i.x, r_i.y, r_i.z, F.x, F.y, F.z);
        }
//        COut::Printf("\r\n\r\n------------------ Minus PBCz -----------------------\r\n");
//        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x, r_k.y, r_k.z, F.x, F.y, F.z, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x - PBCz;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
//            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.Abs() > 10.0) ? "***" : "   ", i, r_i.x, r_i.y, r_i.z, F.x, F.y, F.z);
        }
//        COut::CloseOutputFile();
        exit(0);
    }
}
*/
