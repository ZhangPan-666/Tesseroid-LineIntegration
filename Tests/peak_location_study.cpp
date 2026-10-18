#include "../C++Codes/_TLI_estimate.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Tesseroid::Component;
using Tesseroid::Edge;
using Tesseroid::EdgeKind;
using Tesseroid::Point;
using Tesseroid::Prism;

constexpr double rad(double x) { return x*Tesseroid::PI/180.0; }

struct PeakSummary {
    double location=0.0, amplitude=0.0, signed_value=0.0;
    int peaks_1e2=0, peaks_1e3=0;
};

template<class F>
static double golden_abs(const F& f,double a,double b) {
    constexpr double q=0.6180339887498948482;
    double c=b-q*(b-a),d=a+q*(b-a),yc=std::abs(f(c)),yd=std::abs(f(d));
    for(int i=0;i<80;++i){
        if(yc<yd){a=c;c=d;yc=yd;d=a+q*(b-a);yd=std::abs(f(d));}
        else{b=d;d=c;yd=yc;c=b-q*(b-a);yc=std::abs(f(c));}
    }
    return yc>yd?c:d;
}

static std::vector<double> scan_grid(double a,double b,const std::vector<double>& centres,
                                     double scale) {
    std::vector<double> x;
    x.reserve(8001+centres.size()*3001);
    for(int i=0;i<=8000;++i)x.push_back(a+(b-a)*i/8000.0);
    for(double c:centres){
        x.push_back(c);
        for(int i=0;i<1500;++i){
            const double q=std::pow(10.0,-9.0+12.0*i/1499.0);
            const double d=std::max(scale,1e-15*(b-a))*q;
            if(c-d>a)x.push_back(c-d);
            if(c+d<b)x.push_back(c+d);
        }
    }
    std::sort(x.begin(),x.end());
    x.erase(std::unique(x.begin(),x.end()),x.end());
    return x;
}

template<class F>
static PeakSummary summarize(const F& f,const std::vector<double>& x) {
    std::vector<double> y(x.size()),ay(x.size());
    for(std::size_t i=0;i<x.size();++i){y[i]=f(x[i]);ay[i]=std::abs(y[i]);}
    const auto imax=static_cast<std::size_t>(std::max_element(ay.begin(),ay.end())-ay.begin());
    struct Candidate{double x,a,y;};
    std::vector<Candidate> peaks;
    auto add=[&](std::size_t i,double,double){
        peaks.push_back({x[i],ay[i],y[i]});
    };
    if(ay.front()>=ay[1])peaks.push_back({x.front(),ay.front(),y.front()});
    for(std::size_t i=1;i+1<x.size();++i)
        if(ay[i]>=ay[i-1]&&ay[i]>=ay[i+1])add(i,x[i-1],x[i+1]);
    if(ay.back()>=ay[x.size()-2])peaks.push_back({x.back(),ay.back(),y.back()});
    if(peaks.empty())add(imax,x[imax?imax-1:imax],x[std::min(imax+1,x.size()-1)]);
    std::sort(peaks.begin(),peaks.end(),[](const auto& p,const auto& q){return p.a>q.a;});
    const double top=peaks.front().a;
    PeakSummary out{peaks.front().x,top,peaks.front().y,0,0};
    for(const auto& p:peaks){
        if(p.a>=1e-2*top)++out.peaks_1e2;
        if(p.a>=1e-3*top)++out.peaks_1e3;
    }
    return out;
}

static double angular_distance(const Point& p,double lat,double lon) {
    const double c=std::sin(p.lat)*std::sin(lat)+
        std::cos(p.lat)*std::cos(lat)*std::cos(lon-p.lon);
    return std::acos(std::clamp(c,-1.0,1.0));
}

int main(int argc,char** argv) {
    const fs::path output=argc>1?fs::path(argv[1]):
        fs::path("results/peak_location_study.csv");
    fs::create_directories(output.parent_path());
    std::ofstream out(output);out<<std::setprecision(17);
    out<<"cell_width_deg,kernel_form,edge_kind,component,height_km,gap_deg,along_fraction,"
           "coordinate_projection_rad,geometric_projection_rad,global_peak_rad,"
           "peak_minus_geometric_rad,peak_minus_coordinate_rad,offset_over_core,"
           "peak_amplitude,peak_signed_value,peaks_above_1e2,peaks_above_1e3\n";

    const std::array<Component,4> components{
        Component::V,Component::Vz,Component::Vzz,Component::Vzzz};
    for(double width_deg:{0.1,20.0}){
        const double latc=rad(35.0),lonc=rad(120.0),h=0.5*rad(width_deg);
        const Prism prism{latc-h,latc+h,lonc-h,lonc+h,6341.2,6371.2,1.0};
        for(EdgeKind kind:{EdgeKind::Latitude,EdgeKind::Longitude}){
            const double a=kind==EdgeKind::Latitude?prism.lon1:prism.lat1;
            const double b=kind==EdgeKind::Latitude?prism.lon2:prism.lat2;
            for(double height:{0.001,10.0})for(double gap:{1e-6,1e-3,0.1})
            for(double fraction:{0.1,0.5,0.9}){
                Point point{};point.r=prism.r2+height;
                double fixed=0.0,coordinate=0.0;
                if(kind==EdgeKind::Latitude){
                    fixed=prism.lat2;point.lat=fixed+rad(gap);
                    point.lon=a+fraction*(b-a);coordinate=point.lon;
                }else{
                    fixed=prism.lon2;point.lon=fixed+rad(gap);
                    point.lat=a+fraction*(b-a);coordinate=point.lat;
                }
                const Edge edge{kind,fixed,a,b,+1};
                const double geometric=Tesseroid::detail::geometric_projection(edge,point);
                const double gamma=kind==EdgeKind::Latitude?
                    angular_distance(point,fixed,geometric):angular_distance(point,geometric,fixed);
                const double core=std::hypot(gamma,height/point.r);
                const auto grid=scan_grid(a,b,{coordinate,geometric},core);
                for(Component component:components){
                    auto single=[&](double x){return Tesseroid::kernel(component,prism,point,edge,x);};
                    auto paired=[&](double x){return Tesseroid::paired_kernel(component,prism,point,kind,x);};
                    for(const auto& form:std::array<std::pair<const char*,const void*>,2>{
                            std::pair<const char*,const void*>{"single_edge",&single},
                            std::pair<const char*,const void*>{"paired_edges",&paired}}){
                        PeakSummary s=form.first[0]=='s'?summarize(single,grid):summarize(paired,grid);
                        out<<width_deg<<','<<form.first<<','
                           <<(kind==EdgeKind::Latitude?"latitude":"longitude")<<','
                           <<Tesseroid::component_name(component)<<','<<height<<','<<gap<<','<<fraction<<','
                           <<coordinate<<','<<geometric<<','<<s.location<<','
                           <<s.location-geometric<<','<<s.location-coordinate<<','
                           <<std::abs(s.location-geometric)/std::max(core,1e-300)<<','
                           <<s.amplitude<<','<<s.signed_value<<','<<s.peaks_1e2<<','<<s.peaks_1e3<<'\n';
                    }
                }
            }
        }
    }
    std::cout<<"wrote "<<output.string()<<'\n';
}
