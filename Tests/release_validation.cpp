#include "_TLI_estimate.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
using A = std::array<double,20>;
struct Node { double x, w; };
// Independent tensor-product Gauss-Legendre volume integration of Cartesian
// Newton derivatives. No TLI kernels or TLI quadrature routines are used.
std::vector<Node> gauss(int n, double a, double b) {
    std::vector<Node> nodes(n);
    for (int i=0; i<(n+1)/2; ++i) {
        double z=std::cos(Tesseroid::PI*(i+0.75)/(n+0.5)), dp=0;
        for(int it=0; it<50; ++it) {
            double p=1, prev=0;
            for(int j=1; j<=n; ++j) { double old=prev; prev=p; p=((2*j-1)*z*prev-(j-1)*old)/j; }
            dp=n*(z*p-prev)/(z*z-1);
            double next=z-p/dp;
            if(std::abs(next-z)<2e-16) { z=next; break; }
            z=next;
        }
        const double w=(b-a)/((1-z*z)*dp*dp), mid=(a+b)/2, half=(b-a)/2;
        nodes[i]={mid-half*z,w}; nodes[n-1-i]={mid+half*z,w};
    }
    return nodes;
}
A volume(const Tesseroid::Prism& p, const Tesseroid::Point& x, int n) {
    const auto rn=gauss(n,p.r1,p.r2), ln=gauss(n,p.lon1,p.lon2), fn=gauss(n,p.lat1,p.lat2);
    const std::array<double,3> north{-std::sin(x.lat)*std::cos(x.lon),-std::sin(x.lat)*std::sin(x.lon),std::cos(x.lat)};
    const std::array<double,3> east{-std::sin(x.lon),std::cos(x.lon),0};
    const std::array<double,3> up{std::cos(x.lat)*std::cos(x.lon),std::cos(x.lat)*std::sin(x.lon),std::sin(x.lat)};
    const std::array<std::array<double,3>,3> frame{north,east,up};
    constexpr int ij[6][2]={{0,0},{0,1},{0,2},{1,1},{1,2},{2,2}};
    constexpr int ijk[10][3]={{0,0,0},{0,0,1},{0,0,2},{0,1,1},{0,1,2},{0,2,2},{1,1,1},{1,1,2},{1,2,2},{2,2,2}};
    std::array<long double,20> sum{};
    for(auto r:rn) for(auto f:fn) for(auto l:ln) {
        const std::array<double,3> y{r.x*std::cos(f.x)*std::cos(l.x),r.x*std::cos(f.x)*std::sin(l.x),r.x*std::sin(f.x)};
        std::array<double,3> q{};
        for(int k=0;k<3;++k) for(int j=0;j<3;++j) q[k]+=frame[k][j]*(x.r*up[j]-y[j]);
        const double d2=q[0]*q[0]+q[1]*q[1]+q[2]*q[2], inv=1/std::sqrt(d2);
        const double i3=inv/d2, i5=i3/d2, i7=i5/d2;
        A v{}; v[0]=inv;
        for(int k=0;k<3;++k) v[k+1]=-q[k]*i3;
        for(int k=0;k<6;++k) {int i=ij[k][0],j=ij[k][1]; v[k+4]=3*q[i]*q[j]*i5-(i==j?i3:0);}
        for(int k=0;k<10;++k) {int i=ijk[k][0],j=ijk[k][1],h=ijk[k][2]; v[k+10]=-15*q[i]*q[j]*q[h]*i7+3*i5*((i==j?q[h]:0)+(i==h?q[j]:0)+(j==h?q[i]:0));}
        const long double w=r.w*f.w*l.w*r.x*r.x*std::cos(f.x)*p.density;
        for(int k=0;k<20;++k) sum[k]+=w*v[k];
    }
    A out{}; for(int k=0;k<20;++k) out[k]=static_cast<double>(sum[k]); return out;
}
void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
}
int main() try {
    using namespace Tesseroid;
    constexpr double deg=PI/180;
    const Prism p{31.3*deg,34.7*deg,112.2*deg,116.6*deg,6341.2,6371.2,1.7};
    EstimateOptions o; o.gravitational_constant=1; o.integration.atol=1e-16; o.integration.rtol=1e-10;
    // Nonsymmetric exterior points: independent evidence for every component.
    for(const Point x : {Point{33.11*deg,114.08*deg,6771.2}, Point{28.4*deg,118.2*deg,6471.2}}) {
        const auto low=volume(p,x,20), high=volume(p,x,32);
        const auto got=estimate_all(p,x,o);
        double worst=0;
        for(int k=0;k<20;++k) {
            const double scale=std::max(std::abs(high[k]),1e-25);
            require(std::abs(low[k]-high[k])<2e-8*scale,"volume reference did not converge");
            const double error=std::abs(got[k].value-high[k])/scale;
            std::cout << component_name(static_cast<Component>(k)) << " relative_error=" << error << '\n';
            require(error<2e-7,"TLI disagrees with volume reference");
            require(got[k].converged,"TLI did not meet its requested tolerance");
            worst=std::max(worst,error);
        }
        std::cout << "all-20 maximum relative error: " << worst << '\n';
    }
    const Point interior{33.11*deg,114.08*deg,6357.1};
    const auto inside=estimate_all(p,interior,o);
    require(std::abs(inside[4].value+inside[7].value+inside[9].value+4*PI*p.density)<1e-7,"interior Poisson trace failed");
    Prism zero=p; zero.density=0;
    for(const auto& v:estimate_all(zero,Point{33*deg,114*deg,6771.2},o)) require(v.value==0,"zero density failed");
    bool rejected=false;
    try { Prism invalid=p; invalid.r1=invalid.r2; estimate_all(invalid,interior,o); } catch(const std::invalid_argument&) {rejected=true;}
    require(rejected,"invalid geometry accepted");
    rejected=false;
    try { Point face=interior; face.r=p.r1; estimate_all(p,face,o); } catch(const std::domain_error&) {rejected=true;}
    require(rejected,"undefined lower-face convention accepted");
    std::cout << "Independent volume, Poisson, zero-density and geometry checks passed\n";
    return 0;
} catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
