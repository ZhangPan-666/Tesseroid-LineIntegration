#include "../C++Codes/_TLI_estimate.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace fs=std::filesystem;
using namespace Tesseroid;
constexpr double rad(double x){return x*PI/180.0;}

struct Result{double x=0,a=0;int signed_max=0,signed_min=0,abs_max=0;};

static Result inspect(const std::vector<double>& x,const std::vector<double>& y){
    Result r;std::vector<double> a(y.size());
    std::transform(y.begin(),y.end(),a.begin(),[](double z){return std::abs(z);});
    const auto im=std::max_element(a.begin(),a.end())-a.begin();r.x=x[im];r.a=a[im];
    const double level=.01*r.a,eps=1e-10*r.a;
    for(std::size_t i=1;i+1<y.size();++i){
        if(std::abs(y[i])<level)continue;
        if(y[i]>=y[i-1]&&y[i]>=y[i+1]&&
           (y[i]>y[i-1]+eps||y[i]>y[i+1]+eps))++r.signed_max;
        if(y[i]<=y[i-1]&&y[i]<=y[i+1]&&
           (y[i]<y[i-1]-eps||y[i]<y[i+1]-eps))++r.signed_min;
        if(a[i]>=a[i-1]&&a[i]>=a[i+1]&&
           (a[i]>a[i-1]+eps||a[i]>a[i+1]+eps))++r.abs_max;
    }
    return r;
}

int main(int argc,char**argv){
    const fs::path output=argc>1?fs::path(argv[1]):fs::path("results/radial_peak_profiles.csv");
    fs::create_directories(output.parent_path());std::ofstream out(output);out<<std::setprecision(17);
    out<<"edge_kind,kernel_form,component,height_km,gap_deg,geometric_rad,coordinate_rad,"
           "global_peak_rad,peak_offset_over_core,signed_maxima,signed_minima,abs_maxima\n";
    const double h=.5*rad(.1),latc=rad(35),lonc=rad(120);
    const Prism prism{latc-h,latc+h,lonc-h,lonc+h,6341.2,6371.2,1};
    const std::array<int,4> ci{0,3,9,19};
    for(EdgeKind kind:{EdgeKind::Latitude,EdgeKind::Longitude})
    for(const auto hg:std::array<std::pair<double,double>,4>{{{.001,1e-6},{.001,1e-3},{.1,1e-3},{10,.1}}}){
        const double a=kind==EdgeKind::Latitude?prism.lon1:prism.lat1;
        const double b=kind==EdgeKind::Latitude?prism.lon2:prism.lat2;
        const double fixed=kind==EdgeKind::Latitude?prism.lat2:prism.lon2;
        Point point{};point.r=prism.r2+hg.first;
        if(kind==EdgeKind::Latitude){point.lat=fixed+rad(hg.second);point.lon=.5*(a+b);}
        else{point.lon=fixed+rad(hg.second);point.lat=.5*(a+b);}
        Edge edge{kind,fixed,a,b,+1};
        const double coordinate=kind==EdgeKind::Latitude?point.lon:point.lat;
        const double geometric=detail::geometric_projection(edge,point);
        const double gamma=detail::angular_separation(point,
            kind==EdgeKind::Latitude?fixed:geometric,kind==EdgeKind::Latitude?geometric:fixed);
        const double core=std::hypot(gamma,hg.first/point.r);
        double la=std::max(a,geometric-50*core),lb=std::min(b,geometric+50*core);
        std::vector<double>x; x.reserve(40005);
        for(int i=0;i<=40000;++i)x.push_back(la+(lb-la)*i/40000.0);
        x.push_back(a);x.push_back(b);x.push_back(geometric);x.push_back(coordinate);
        std::sort(x.begin(),x.end());x.erase(std::unique(x.begin(),x.end()),x.end());
        std::array<std::vector<double>,20> single,paired;
        for(auto&v:single)v.reserve(x.size());for(auto&v:paired)v.reserve(x.size());
        for(double z:x){
            const auto s=kernel_all(prism,point,edge,z),p=paired_kernel_all(prism,point,kind,z);
            for(int k:ci){single[k].push_back(s[k]);paired[k].push_back(p[k]);}
        }
        for(int k:ci)for(int form=0;form<2;++form){
            const auto r=inspect(x,form?paired[k]:single[k]);
            out<<(kind==EdgeKind::Latitude?"latitude":"longitude")<<','
               <<(form?"paired_edges":"single_edge")<<','
               <<component_name(static_cast<Component>(k))<<','<<hg.first<<','<<hg.second<<','
               <<geometric<<','<<coordinate<<','<<r.x<<','
               <<std::abs(r.x-geometric)/std::max(core,1e-300)<<','
               <<r.signed_max<<','<<r.signed_min<<','<<r.abs_max<<'\n';
        }
    }
    std::cout<<"wrote "<<output.string()<<'\n';
}
