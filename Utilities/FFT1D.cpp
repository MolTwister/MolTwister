#include <math.h>
#include <memory>
#include "FFT1D.h"

std::shared_ptr<std::vector<CFFT1D::CCplx>> CFFT1D::fft1D(const std::vector<CCplx> &in, EDir dir) const
{
    std::shared_ptr<std::vector<CCplx>> omega = buildOmega((int)in.size(), dir);

    auto out = std::make_shared<std::vector<CCplx>>(in.size());
    fft1DRec(in, 0, *out, (int)in.size(), *omega, 1, dir);

    return out;
}

std::shared_ptr<std::vector<CFFT1D::CCplx>> CFFT1D::dft1D(const std::vector<CCplx> &in, int stride, EDir) const
{
    auto out = std::make_shared<std::vector<CCplx>>(in.size());
    dft1D(in, 0, static_cast<int>(in.size()), *out, stride);

    return out;
}

std::shared_ptr<std::vector<CFFT1D::CCplx>> CFFT1D::buildOmega(int N, EDir dir) const
{
    double dirNum = double(dir);
    
    auto omega = std::make_shared<std::vector<CCplx>>(N*N);
    
    for(int i=0; i<N; i++)
    {
        for(int n=0; n<N; n++)
        {
            (*omega)[n + i*N].re_ = cos(dirNum*2.0*M_PI * double(i) / double(n+1));
            (*omega)[n + i*N].im_ = sin(dirNum*2.0*M_PI * double(i) / double(n+1));
        }
    }

    return omega;
}

void CFFT1D::dft1D(const std::vector<CCplx> &in, int offset, int n, std::vector<CCplx>& out, int stride, EDir dir) const
{
    float u;
    float dirNum = float(dir);
    
    out.clear();
    out.resize(n);
    
    for(int k=0; k<n; k++)
    {
        for(int j=0; j<n; j++)
        {
            u = dirNum*2.0f*float(M_PI)*float(j)*float(k)/float(n);
            out[k].re_+= in[offset + j*stride].re_*cos((double)u) - in[offset + j*stride].im_*sin((double)u);
            out[k].im_+= in[offset + j*stride].re_*sin((double)u) + in[offset + j*stride].im_*cos((double)u);
        }
    }
}

void CFFT1D::fft1DRec(const std::vector<CCplx> &in, int offset, std::vector<CCplx>& out, int n, std::vector<CCplx>& omega, int stride, EDir dir) const
{
    std::vector<CCplx> ye, yo;
    double tauRe, tauIm;
    int m;

    ye.clear();
    ye.resize(n/2);
    yo.clear();
    yo.resize(n/2);
    
    if(n % 2 == 0)
    {
        fft1DRec(in, offset, ye, n/2, omega, stride*2, dir);
        fft1DRec(in, offset+stride, yo, n/2, omega, stride*2, dir);
        
        for(int k=0; k<n/2; k++)
        {
            m = n + int(in.size())*k - 1;
            tauRe = omega[m].re_*yo[k].re_ - omega[m].im_*yo[k].im_;
            tauIm = omega[m].re_*yo[k].im_ + omega[m].im_*yo[k].re_;
            
            out[k].re_ = ye[k].re_ + tauRe;
            out[k].im_ = ye[k].im_ + tauIm;
            out[k+n/2].re_ = ye[k].re_ - tauRe;
            out[k+n/2].im_ = ye[k].im_ - tauIm;
        }
    }
    else dft1D(in, offset, n, out, stride, dir);
}

void CFFT1D::zeroPad(std::vector<CCplx>& data)
{
    int newNum = 2;

    while(newNum < data.size()) newNum*= 2;

    for(int i=(int)data.size(); i<newNum; i++)
    {
        CCplx zero;
        data.emplace_back(zero);
    }
}
