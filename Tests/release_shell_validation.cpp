#include "_tesseroid_MSI.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
    using namespace Tesseroid;
    constexpr double a=6341.2,b=6371.2;
    const auto grid=tli_mmsi_detail::make_grid(-175,175,-85,85,10,10,648);
    const std::vector<double> outer(648,b),inner(648,a),rho(648,1);
    for(double r : {6771.2,6300.0,6356.2}) {
        const auto value=tli_mmsi_detail::point_estimate(grid,outer,inner,rho,24.37,138.91,r,1e-12,1e-9);
        double V,dV,ddV,dddV;
        if(r>=b) {
            const double mass=4*PI*(b*b*b-a*a*a)/3;
            V=mass/r; dV=-mass/(r*r); ddV=2*mass/(r*r*r); dddV=-6*mass/std::pow(r,4);
        } else if(r<=a) {V=2*PI*(b*b-a*a);dV=ddV=dddV=0;}
        else {
            V=4*PI*((r*r*r-a*a*a)/(3*r)+(b*b-r*r)/2);
            dV=-4*PI*(r-a*a*a/(r*r))/3;
            ddV=-4*PI*(1+2*a*a*a/(r*r*r))/3;
            dddV=8*PI*a*a*a/std::pow(r,4);
        }
        std::array<double,20> ref{};
        ref[0]=V;ref[3]=dV;ref[4]=ref[7]=dV/r;ref[9]=ddV;
        ref[12]=ref[17]=ddV/r-dV/(r*r);ref[19]=dddV;
        const double order_scale[4]={std::max(std::abs(V),1.0),std::max(std::abs(dV),1.0),std::max(std::abs(ddV),1.0),std::max(std::abs(dddV),1e-3)};
        for(int k=0;k<20;++k) {
            const double err=std::abs(value[k]-ref[k]);
            const int order=component_order(static_cast<Component>(k));
            if(!std::isfinite(value[k]) || err>2e-7*order_scale[order]) {
                std::cerr << "shell r=" << r << ' ' << component_name(static_cast<Component>(k)) << " error=" << err << '\n';return 1;
            }
        }
        std::cout << "Closed-form shell reference passed at radius " << r << '\n';
    }
    return 0;
}
