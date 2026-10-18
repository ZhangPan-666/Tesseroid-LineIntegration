#include "../C++Codes/_TLI_estimate.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Tesseroid::Component;
using Tesseroid::EdgeKind;
using Tesseroid::Point;
using Tesseroid::Prism;

constexpr double rad(double degree) { return degree*Tesseroid::PI/180.0; }

struct Geometry {
    Point point;
    double a=0,b=0,t0=0,gamma=0,distance_km=0,core=0;
};

struct Features {
    double location=0,amplitude=0,fwhm=0,support=0,asymmetry=0,l1=0;
};

struct MethodResult {
    double value=0,error=std::numeric_limits<double>::quiet_NaN(),milliseconds=0;
    std::size_t evaluations=0,intervals=0;
    bool converged=true;
};

static Geometry geometry(const Prism& p,EdgeKind kind,double height_km,double gap_deg,
                         double along_fraction) {
    Geometry g;
    g.point.r=p.r2+height_km;
    if(kind==EdgeKind::Latitude){
        g.a=p.lon1; g.b=p.lon2; g.t0=g.a+along_fraction*(g.b-g.a);
        g.point.lat=p.lat2+rad(gap_deg); g.point.lon=g.t0;
        g.gamma=std::abs(g.point.lat-p.lat2);
    }else{
        g.a=p.lat1; g.b=p.lat2; g.point.lat=g.a+along_fraction*(g.b-g.a);
        g.point.lon=p.lon2+rad(gap_deg);
        g.t0=std::atan2(std::sin(g.point.lat),
            std::cos(g.point.lat)*std::cos(p.lon2-g.point.lon));
        g.t0=std::clamp(g.t0,g.a,g.b);
        const double c=std::sin(g.point.lat)*std::sin(g.t0)+
            std::cos(g.point.lat)*std::cos(g.t0)*std::cos(g.point.lon-p.lon2);
        g.gamma=std::acos(std::clamp(c,-1.0,1.0));
    }
    const double radial=std::max(0.0,g.point.r-p.r2);
    g.distance_km=std::sqrt(radial*radial+2*g.point.r*p.r2*(1-std::cos(g.gamma)));
    g.core=std::hypot(radial/g.point.r,g.gamma);
    return g;
}

static std::pair<std::vector<double>,std::vector<double>> gauss_legendre(int n) {
    std::vector<double> x(n),w(n);
    const int m=(n+1)/2;
    for(int i=0;i<m;++i){
        double z=std::cos(Tesseroid::PI*(i+0.75)/(n+0.5)),z1;
        double pp=0;
        do{
            double p1=1,p2=0;
            for(int j=1;j<=n;++j){
                const double p3=p2;p2=p1;
                p1=((2*j-1)*z*p2-(j-1)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1);
            z1=z;z=z1-p1/pp;
        }while(std::abs(z-z1)>2e-15);
        x[i]=-z;x[n-1-i]=z;
        w[i]=2/((1-z*z)*pp*pp);w[n-1-i]=w[i];
    }
    return {x,w};
}

template<class F>
static std::array<double,20> gauss_vector(const F& f,double a,double b,int n) {
    const auto [x,w]=gauss_legendre(n);
    std::array<double,20> out{};
    const double c=0.5*(a+b),h=0.5*(b-a);
    for(int i=0;i<n;++i){
        const auto y=f(c+h*x[i]);
        for(int k=0;k<20;++k)out[k]=std::fma(h*w[i],y[k],out[k]);
    }
    return out;
}

template<class F>
static std::array<double,20> composite_gauss(const F& f,const std::vector<double>& cuts,int n) {
    std::array<double,20> out{};
    for(std::size_t i=1;i<cuts.size();++i){
        const auto q=gauss_vector(f,cuts[i-1],cuts[i],n);
        for(int k=0;k<20;++k)out[k]+=q[k];
    }
    return out;
}

static std::vector<double> sampling_grid(double a,double b,double t0,double core,
                                         int uniform_count,int local_count) {
    std::vector<double> x;
    x.reserve(uniform_count+2*local_count+3);
    for(int i=0;i<uniform_count;++i)x.push_back(a+(b-a)*i/(uniform_count-1.0));
    const double width=std::max(core,1e-14*(b-a));
    for(int i=0;i<local_count;++i){
        const double q=std::pow(10.0,-7.0+10.0*i/std::max(1,local_count-1));
        x.push_back(t0-width*q);x.push_back(t0+width*q);
    }
    x.push_back(t0);
    x.erase(std::remove_if(x.begin(),x.end(),[&](double z){return z<a||z>b;}),x.end());
    std::sort(x.begin(),x.end());
    x.erase(std::unique(x.begin(),x.end()),x.end());
    return x;
}

static double crossing(const std::vector<double>& x,const std::vector<double>& y,
                       std::size_t i0,int direction,double level) {
    std::size_t i=i0;
    while((direction<0&&i>0)||(direction>0&&i+1<x.size())){
        const std::size_t j=static_cast<std::size_t>(static_cast<long long>(i)+direction);
        if(y[j]<=level){
            const double den=y[i]-y[j];
            return den==0?x[j]:x[j]+(level-y[j])*(x[i]-x[j])/den;
        }
        i=j;
    }
    return direction<0?x.front():x.back();
}

static Features features(const std::vector<double>& x,const std::vector<double>& y) {
    Features f;
    std::vector<double> ay(y.size());
    std::transform(y.begin(),y.end(),ay.begin(),[](double q){return std::abs(q);});
    const auto peak=static_cast<std::size_t>(std::max_element(ay.begin(),ay.end())-ay.begin());
    f.location=x[peak];f.amplitude=ay[peak];
    if(f.amplitude==0)return f;
    const double left=crossing(x,ay,peak,-1,0.5*f.amplitude);
    const double right=crossing(x,ay,peak,1,0.5*f.amplitude);
    f.fwhm=right-left;
    const double lw=f.location-left,rw=right-f.location;
    f.asymmetry=(rw-lw)/std::max(rw+lw,std::numeric_limits<double>::min());
    const double support_left=crossing(x,ay,peak,-1,1e-3*f.amplitude);
    const double support_right=crossing(x,ay,peak,1,1e-3*f.amplitude);
    f.support=support_right-support_left;
    for(std::size_t i=1;i<x.size();++i)
        f.l1+=0.5*(x[i]-x[i-1])*(ay[i]+ay[i-1]);
    return f;
}

static std::vector<double> reference_cuts(double a,double b,double t0,double core) {
    auto cuts=sampling_grid(a,b,t0,core,129,100);
    return cuts;
}

template<class F>
static MethodResult adaptive(const F& f,double a,double b,double atol,double rtol,
                             bool use_gk21=false) {
    Tesseroid::IntegralOptions o;
    o.atol=atol;o.rtol=rtol;o.max_evaluations=100000;o.max_intervals=5000;
    o.use_gk21=use_gk21;
    Tesseroid::IntegralInfo info;
    const auto start=Clock::now();
    const double value=Tesseroid::integral(f,a,b,o,&info);
    const auto stop=Clock::now();
    return {value,info.error,std::chrono::duration<double,std::milli>(stop-start).count(),
        info.evaluations,info.intervals,info.converged};
}

template<class F>
static MethodResult split_adaptive(const F& f,double a,double b,const Features& feature,
                                   double core,double atol,double rtol) {
    std::vector<double> cuts{a,b,feature.location};
    const double hw=std::max(0.5*feature.fwhm,core);
    for(double m:{0.25,0.5,1.0,2.0,4.0,8.0}){
        cuts.push_back(feature.location-m*hw);cuts.push_back(feature.location+m*hw);
    }
    cuts.erase(std::remove_if(cuts.begin(),cuts.end(),[&](double z){return z<=a||z>=b;}),cuts.end());
    cuts.push_back(a);cuts.push_back(b);std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    MethodResult total;total.error=0;total.converged=true;
    const auto start=Clock::now();
    for(std::size_t i=1;i<cuts.size();++i){
        auto q=adaptive(f,cuts[i-1],cuts[i],atol,rtol,true);
        total.value+=q.value;total.error+=q.error;total.evaluations+=q.evaluations;
        total.intervals+=q.intervals;total.converged=total.converged&&q.converged;
    }
    total.milliseconds=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
    return total;
}

static void write_row(std::ofstream& out,int case_id,const char* kind,Component component,
    double height,double gap,double along_fraction,const Geometry& g,const Features& f,const char* method,
    const MethodResult& q,double reference,double reference_error) {
    const double actual=std::abs(q.value-reference);
    const double rel=actual/std::max(std::abs(reference),std::numeric_limits<double>::min());
    const double ratio=std::isfinite(q.error)?actual/std::max(q.error,std::numeric_limits<double>::min()):
        std::numeric_limits<double>::quiet_NaN();
    out<<case_id<<','<<kind<<','<<Tesseroid::component_name(component)<<','
       <<Tesseroid::component_order(component)<<','<<height<<','<<gap<<','<<along_fraction<<','
       <<g.gamma<<','<<g.distance_km<<','<<g.core<<','<<g.core/(g.b-g.a)<<','
       <<f.location<<','<<f.amplitude<<','<<f.fwhm<<','<<f.fwhm/(g.b-g.a)<<','
       <<f.support<<','<<f.support/(g.b-g.a)<<','<<f.asymmetry<<','<<f.l1<<','
       <<method<<','<<q.value<<','<<q.error<<','<<actual<<','<<rel<<','<<ratio<<','
       <<q.evaluations<<','<<q.intervals<<','<<(q.converged?1:0)<<','<<q.milliseconds<<','
       <<reference<<','<<reference_error<<'\n';
}

int main(int argc,char** argv) {
    const std::string mode=argc>1?argv[1]:"smoke";
    fs::path output=argc>2?fs::path(argv[2]):fs::path("results/kernel_peak_regimes_cpp.csv");
    fs::create_directories(output.parent_path());
    std::ofstream out(output);out<<std::setprecision(17);
    out<<"case_id,kind,component,derivative_order,height_km,angular_gap_deg,along_fraction,minimum_angle_rad,"
          "minimum_distance_km,normalized_distance,predicted_width_ratio,peak_location_rad,"
          "peak_amplitude,fwhm_rad,fwhm_ratio,support_width_rad,support_width_ratio,asymmetry,"
          "l1_area,method,value,estimated_error,actual_abs_error,actual_rel_error,estimator_ratio,"
          "function_evaluations,intervals,converged,time_ms,reference_value,reference_error\n";

    const Prism prism{rad(30),rad(50),rad(110),rad(130),6341.2,6371.2,1.0};
    std::vector<int> components{0,3,9,19};
    std::vector<double> heights{0.01,15},gaps{1e-4,1e-2,1};
    std::vector<double> along_fractions{0.95};
    int uniform_samples=2001,local_samples=300;
    if(mode=="quick"){
        components={0,3,4,9,10,14,19};heights={0.01,1,15,100};
        gaps={1e-4,5e-4,0.0025,0.0125,0.0625,0.3125,1.5};
        uniform_samples=10001;local_samples=1000;
    }else if(mode=="position"){
        components={19};heights={0.01};gaps={1e-4};
        along_fractions={0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
                         0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95};
        uniform_samples=1001;local_samples=200;
    }else if(mode=="full"){
        components.resize(20);std::iota(components.begin(),components.end(),0);
        heights={0.001,0.01,0.1,1,15,100,400};
        gaps={1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.3,1,3,10};
        uniform_samples=20001;local_samples=2000;
    }else if(mode!="smoke"){
        std::cerr<<"mode must be smoke, quick, full, or position\n";return 2;
    }

    int case_id=0;
    for(EdgeKind kind:{EdgeKind::Latitude,EdgeKind::Longitude}){
        const char* kind_name=kind==EdgeKind::Latitude?"latitude":"longitude";
        for(double height:heights)for(double gap:gaps)for(double along_fraction:along_fractions){
            const Geometry g=geometry(prism,kind,height,gap,along_fraction);
            Tesseroid::detail::PreparedPairedKernel prepared;
            { const Tesseroid::detail::ScopedJetContext context(20);
              prepared=Tesseroid::prepare_paired_kernel(prism,g.point,kind); }
            auto all=[&](double x){
                const Tesseroid::detail::ScopedJetContext context(20);
                return Tesseroid::detail::derivatives(
                    Tesseroid::paired_kernel_jet(prism,prepared,x));
            };

            const auto scan=sampling_grid(g.a,g.b,g.t0,g.core,uniform_samples,local_samples);
            std::array<std::vector<double>,20> sampled;
            for(auto& y:sampled)y.reserve(scan.size());
            for(double x:scan){const auto y=all(x);for(int k=0;k<20;++k)sampled[k].push_back(y[k]);}
            const auto cuts=reference_cuts(g.a,g.b,g.t0,g.core);
            const auto ref32=composite_gauss(all,cuts,32);
            const auto ref24=composite_gauss(all,cuts,24);
            std::array<std::array<double,20>,4> fixed{};
            const int orders[4]={6,10,20,30};
            for(int i=0;i<4;++i)fixed[i]=gauss_vector(all,g.a,g.b,orders[i]);

            for(int ci:components){
                const auto component=static_cast<Component>(ci);++case_id;
                const Features feature=features(scan,sampled[ci]);
                const double reference=ref32[ci],reference_error=std::abs(ref32[ci]-ref24[ci]);
                const double scale=std::max(feature.l1,feature.amplitude*std::max(feature.fwhm,1e-30));
                const double atol=std::max(1e-10*scale,std::numeric_limits<double>::min());
                MethodResult ref{reference,reference_error,0,cuts.size()*32,cuts.size()-1,true};
                write_row(out,case_id,kind_name,component,height,gap,along_fraction,g,feature,
                    "ReferenceCompositeGL",ref,reference,reference_error);
                for(int i=0;i<4;++i){
                    MethodResult q{fixed[i][ci],std::numeric_limits<double>::quiet_NaN(),0,
                        static_cast<std::size_t>(orders[i]),1,true};
                    const std::string name="GL"+std::to_string(orders[i]);
                    write_row(out,case_id,kind_name,component,height,gap,along_fraction,g,feature,
                        name.c_str(),q,reference,reference_error);
                }
                auto scalar=[&](double x){return Tesseroid::paired_kernel_prepared(
                    component,prism,prepared,x);};
                const auto ordinary=adaptive(scalar,g.a,g.b,atol,1e-8);
                write_row(out,case_id,kind_name,component,height,gap,along_fraction,g,feature,
                    "AdaptiveGK",ordinary,reference,reference_error);
                if(mode!="position"){
                    const auto split=split_adaptive(scalar,g.a,g.b,feature,g.core,atol,1e-8);
                    write_row(out,case_id,kind_name,component,height,gap,along_fraction,g,feature,
                        "PeakSplitGK",split,reference,reference_error);
                }
                std::cout<<case_id<<' '<<kind_name<<' '<<Tesseroid::component_name(component)
                    <<" h="<<height<<" gap="<<gap<<" FWHM/L="<<feature.fwhm/(g.b-g.a)<<'\n';
            }
        }
    }
    std::cout<<"wrote "<<output.string()<<"\n";
}
