#pragma once
#include <stdio.h>
#include <vector>
#include <memory>

class CFFT1D
{
public:
    enum EDir { dirFwd=-1, dirRev=1 };

public:
    class CCplx
    {
    public:
        CCplx() { re_ = 0.0; im_ = 0.0; }
        CCplx(double re, double im) { re_ = re; im_ = im; }

    public:
        CCplx conj() { return CCplx(re_, -im_); }
        double modulus2() { return re_*re_ + im_*im_; }
        
    public:
        double re_;
        double im_;
    };
    
public:
    CFFT1D() {}
    
public:
    std::shared_ptr<std::vector<CCplx>> fft1D(const std::vector<CCplx>& in, EDir dir=dirFwd) const;
    std::shared_ptr<std::vector<CCplx>> dft1D(const std::vector<CCplx>& in, int stride, EDir dir=dirFwd) const;
    static void zeroPad(std::vector<CCplx>& data);

protected:
    std::shared_ptr<std::vector<CCplx>> buildOmega(int N, EDir dir) const;
    void dft1D(const std::vector<CCplx>& in, int offset, int n, std::vector<CCplx>& out, int stride, EDir dir=dirFwd) const;
    
private:
    void fft1DRec(const std::vector<CCplx>& in, int offset, std::vector<CCplx>& out, int n, std::vector<CCplx>& omega, int stride, EDir dir) const;
};
