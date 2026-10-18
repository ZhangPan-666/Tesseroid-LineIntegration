#pragma once

#include "_integral.h"
#include "_TLI_integralkernel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef TLI_ENABLE_MPI
#include <mpi.h>
#endif

namespace Tesseroid {

struct EstimateOptions {
    IntegralOptions integration{};
    double radial_limit_relative_step = 1.0e-5;
    std::size_t radial_limit_max_refinements = 4;
    double gravitational_constant = 6.67430e-11;
};

struct EstimateInfo {
    Component component = Component::V;
    double value = 0.0;
    double error = 0.0;
    double absolute = 0.0;
    std::size_t integral_evaluations = 0;
    std::size_t intervals = 0;
    bool converged = true;
    double milliseconds = 0.0;
};

namespace detail {

inline double clamp(double x, double a, double b) {
    return std::min(b, std::max(a, x));
}

inline bool close_coordinate(double x,double y){
    return std::abs(x-y)<=64.0*std::numeric_limits<double>::epsilon()*
        std::max({1.0,std::abs(x),std::abs(y)});
}

inline double coordinate_tolerance(double x,double a,double b){
    return 64.0*std::numeric_limits<double>::epsilon()*
        std::max({1.0,std::abs(x),std::abs(a),std::abs(b)});
}

inline bool open_coordinate(double x,double a,double b){
    const double tol=coordinate_tolerance(x,a,b);
    return x>a+tol&&x<b-tol;
}

inline bool closed_coordinate(double x,double a,double b){
    const double tol=coordinate_tolerance(x,a,b);
    return x>=a-tol&&x<=b+tol;
}

enum class BoundaryLocation { Regular, OuterRadialFace, LowerRadialFace, EdgeOrCorner };

inline BoundaryLocation boundary_location(const Prism&p,const Point&x){
    const bool on_r1=close_coordinate(x.r,p.r1);
    const bool on_r2=close_coordinate(x.r,p.r2);
    const bool lat_open=open_coordinate(x.lat,p.lat1,p.lat2);
    const double lon=x.lon+std::round((0.5*(p.lon1+p.lon2)-x.lon)/(2.0*PI))*(2.0*PI);
    const bool lon_open=open_coordinate(lon,p.lon1,p.lon2);
    const bool lat_closed=closed_coordinate(x.lat,p.lat1,p.lat2);
    const bool lon_closed=closed_coordinate(lon,p.lon1,p.lon2);
    if((on_r1||on_r2)&&lat_open&&lon_open)
        return on_r2?BoundaryLocation::OuterRadialFace:BoundaryLocation::LowerRadialFace;
    if((on_r1||on_r2)&&lat_closed&&lon_closed)
        return BoundaryLocation::EdgeOrCorner;

    const bool radial=p.r1<=x.r&&x.r<=p.r2;
    const bool polar=std::abs(std::cos(x.lat))<=64.0*std::numeric_limits<double>::epsilon();
    const bool latitude_face=(close_coordinate(x.lat,p.lat1)||close_coordinate(x.lat,p.lat2))&&
        radial&&(lon_closed||polar);
    const bool longitude_face=(close_coordinate(lon,p.lon1)||close_coordinate(lon,p.lon2))&&
        radial&&lat_closed&&!polar;
    return latitude_face||longitude_face?BoundaryLocation::EdgeOrCorner:BoundaryLocation::Regular;
}

inline double normalize_longitude(double longitude){
    double value=std::remainder(longitude,2.0*PI);
    if(value>=PI)value=-PI;
    return value;
}

inline void validate_geometry(const Prism& prism,const Point& point){
    const bool finite=std::isfinite(prism.r1)&&std::isfinite(prism.r2)&&
        std::isfinite(prism.lat1)&&std::isfinite(prism.lat2)&&
        std::isfinite(prism.lon1)&&std::isfinite(prism.lon2)&&
        std::isfinite(prism.density)&&std::isfinite(point.r)&&
        std::isfinite(point.lat)&&std::isfinite(point.lon);
    if(!finite||!(prism.r1>0.0&&prism.r2>prism.r1&&point.r>0.0&&
       prism.lat1>=-0.5*PI&&prism.lat1<prism.lat2&&prism.lat2<=0.5*PI&&
       prism.lon1>=-PI&&prism.lon1<prism.lon2&&prism.lon2<=PI&&
       point.lat>=-0.5*PI&&point.lat<=0.5*PI&&point.lon>=-PI&&point.lon<PI))
        throw std::invalid_argument(
            "invalid geometry: require positive ordered radii, principal latitude bounds, "
            "and -pi <= lon1 < lon2 <= pi; split antimeridian-crossing tesseroids first");
}

inline double periodic_candidate(double x, double a, double b,
                                      double period) {
    const double k = std::round((0.5*(a+b)-x)/period);
    double best = x + k*period;
    if (best < a) best += period;
    if (best > b) best -= period;
    return clamp(best,a,b);
}

inline double angular_separation(const Point& point,double source_lat,double source_lon){
    const double half_lat=0.5*(source_lat-point.lat);
    const double half_lon=0.5*std::remainder(source_lon-point.lon,2.0*PI);
    const double slat=std::sin(half_lat),slon=std::sin(half_lon);
    const double hav=clamp(slat*slat+std::cos(point.lat)*std::cos(source_lat)*slon*slon,0.0,1.0);
    return 2.0*std::atan2(std::sqrt(hav),std::sqrt(std::max(0.0,1.0-hav)));
}

// Parameter of closest angular approach on one boundary edge.
inline double geometric_projection(const Edge& edge,const Point& point){
    if(edge.kind==EdgeKind::Latitude){
        if(std::abs(std::cos(point.lat))<=
           64.0*std::numeric_limits<double>::epsilon())
            return 0.5*(edge.a+edge.b);
        return periodic_candidate(point.lon,edge.a,edge.b,2.0*PI);
    }
    const double latitude=std::atan2(
        std::sin(point.lat),
        std::cos(point.lat)*std::cos(edge.fixed-point.lon));
    return clamp(latitude,edge.a,edge.b);
}

// At an exact pole, all longitude edges meet at one Cartesian point.  The
// line-integral representation is analytically valid there, but its separate
// edge contributions can be O(gap^-3) while their physical sum is finite.
// Evaluate that removable degeneracy directly in the observation's local
// Cartesian frame.  Longitude is integrated analytically and the remaining
// radial/polar rectangle uses a double-exponential product rule after the
// exact corner-resolving map
//
//   b = 2*sqrt(r*rho)*sin(theta/2),  b = |r-rho|*sinh(s).
//
// Then d^2=(r-rho)^2+b^2 and rho^2*sin(theta)dtheta drho =
// rho*b/r db drho, so even a sub-nanometre radial gap is resolved without
// subtracting large longitude-edge kernels.
struct PolarDENode { double x=0.0,w=0.0; };

inline std::vector<PolarDENode> polar_de_nodes(double step){
    std::vector<PolarDENode> nodes;
    for(int k=0;k<1000;++k){
        const double t=static_cast<double>(k)*step;
        const double u=0.5*PI*std::sinh(t);
        if(u>20.0)break;
        const double ch=std::cosh(u);
        const double x=0.5*(1.0+std::tanh(u));
        const double w=0.25*PI*std::cosh(t)/(ch*ch)*step;
        if(k==0)nodes.push_back({x,w});
        else{
            if(x<1.0)nodes.push_back({x,w});
            const double mirror=1.0-x;
            if(mirror>0.0)nodes.push_back({mirror,w});
        }
    }
    return nodes;
}

// [1,cos,sin,cos^2,cos*sin,sin^2,cos^3,cos^2*sin,cos*sin^2,sin^3]
inline std::array<double,10> polar_longitude_moments(double a,double b){
    const double width=b-a;
    return {{
        width,
        std::sin(b)-std::sin(a),
        std::cos(a)-std::cos(b),
        0.5*width+0.25*(std::sin(2.0*b)-std::sin(2.0*a)),
        0.25*(std::cos(2.0*a)-std::cos(2.0*b)),
        0.5*width-0.25*(std::sin(2.0*b)-std::sin(2.0*a)),
        0.25*(3.0*(std::sin(b)-std::sin(a))+
              (std::sin(3.0*b)-std::sin(3.0*a))/3.0),
        (std::pow(std::cos(a),3)-std::pow(std::cos(b),3))/3.0,
        (std::pow(std::sin(b),3)-std::pow(std::sin(a),3))/3.0,
        (-std::cos(b)+std::pow(std::cos(b),3)/3.0)-
        (-std::cos(a)+std::pow(std::cos(a),3)/3.0)
    }};
}

inline bool polar_cartesian_applicable(const Prism& prism,const Point& point){
    const bool north=close_coordinate(point.lat,0.5*PI);
    const bool south=close_coordinate(point.lat,-0.5*PI);
    const bool radially_separate=point.r>prism.r2||point.r<prism.r1;
    // The Cartesian polar integral is nonsingular whenever the observation
    // and source are separated either radially or angularly.  Requiring only
    // radial separation incorrectly sent exact-r1/r2 observations outside a
    // polar cap through the cancellation-prone line-edge representation.
    const bool angularly_separate=(north||south)&&
        !closed_coordinate(point.lat,prism.lat1,prism.lat2);
    return (north||south)&&(radially_separate||angularly_separate);
}

inline Components polar_cartesian_level(const Prism& prism,const Point& point,
                                        double step,Components* absolute,
                                        std::size_t* evaluations){
    const bool north=point.lat>=0.0;
    const double pole_sign=north?1.0:-1.0;
    const double theta_min=north?0.5*PI-prism.lat2:prism.lat1+0.5*PI;
    const double theta_max=north?0.5*PI-prism.lat1:prism.lat2+0.5*PI;
    const auto moments=polar_longitude_moments(prism.lon1,prism.lon2);
    const auto nodes=polar_de_nodes(step);
    Components total{},correction{},magnitude{};

    // If the observation radius lies inside the source interval, the narrow
    // radial ridge is centred at rho=point.r.  Split there so the DE rule
    // clusters nodes on both sides of the ridge instead of only at r1/r2.
    std::array<std::pair<double,double>,2> radial_panels{{
        {prism.r1,prism.r2},{0.0,0.0}}};
    std::size_t radial_panel_count=1;
    if(point.r>prism.r1&&point.r<prism.r2){
        radial_panels[0]={prism.r1,point.r};
        radial_panels[1]={point.r,prism.r2};
        radial_panel_count=2;
    }

    for(std::size_t panel_index=0;panel_index<radial_panel_count;++panel_index){
      const double radial_a=radial_panels[panel_index].first;
      const double radial_b=radial_panels[panel_index].second;
      const double radial_width=radial_b-radial_a;
      for(const auto& radial_node:nodes)for(const auto& angular_node:nodes){
        const double rho=radial_a+radial_width*radial_node.x;
        const double radial_delta=point.r-rho;
        const double radial_gap=std::abs(radial_delta);
        const double root=std::sqrt(point.r*rho);
        const double b_min=2.0*root*std::sin(0.5*theta_min);
        const double b_max=2.0*root*std::sin(0.5*theta_max);
        double b=0.0;
        double angular_jacobian=0.0;
        if(radial_gap==0.0){
            // A DE radial node can round exactly onto r1/r2 in kilometre
            // coordinates.  When the source is angularly separated from the
            // pole, b_min is positive and direct integration in b is smooth.
            b=b_min+(b_max-b_min)*angular_node.x;
            angular_jacobian=(b_max-b_min)*angular_node.w;
        }else{
            const double s_min=std::asinh(b_min/radial_gap);
            const double s_max=std::asinh(b_max/radial_gap);
            const double source_parameter=
                s_min+(s_max-s_min)*angular_node.x;
            b=radial_gap*std::sinh(source_parameter);
            angular_jacobian=radial_gap*std::cosh(source_parameter)*
                (s_max-s_min)*angular_node.w;
        }
        const double half_sine=b/(2.0*root);
        const double half_cosine=std::sqrt(
            std::max(0.0,1.0-half_sine*half_sine));
        const double transverse=b*std::sqrt(rho/point.r)*half_cosine;
        const double qx_scale=pole_sign*transverse;
        const double qy_scale=-transverse;
        const double qz=radial_delta+b*b/(2.0*point.r);
        const double distance2=radial_delta*radial_delta+b*b;
        const double inverse=1.0/std::sqrt(distance2);
        const double inverse3=inverse/distance2;
        const double inverse5=inverse3/distance2;
        const double inverse7=inverse5/distance2;
        const double jacobian=(rho*b/point.r)*
            radial_width*radial_node.w*
            angular_jacobian;

        Components value{};
        value[0]=moments[0]*inverse;
        value[1]=-qx_scale*moments[1]*inverse3;
        value[2]=-qy_scale*moments[2]*inverse3;
        value[3]=-qz*moments[0]*inverse3;
        value[4]=3.0*qx_scale*qx_scale*moments[3]*inverse5-
                 moments[0]*inverse3;
        value[5]=3.0*qx_scale*qy_scale*moments[4]*inverse5;
        value[6]=3.0*qx_scale*qz*moments[1]*inverse5;
        value[7]=3.0*qy_scale*qy_scale*moments[5]*inverse5-
                 moments[0]*inverse3;
        value[8]=3.0*qy_scale*qz*moments[2]*inverse5;
        value[9]=(3.0*qz*qz*inverse5-inverse3)*moments[0];
        value[10]=9.0*qx_scale*moments[1]*inverse5-
                  15.0*qx_scale*qx_scale*qx_scale*moments[6]*inverse7;
        value[11]=3.0*qy_scale*moments[2]*inverse5-
                  15.0*qx_scale*qx_scale*qy_scale*moments[7]*inverse7;
        value[12]=3.0*qz*moments[0]*inverse5-
                  15.0*qx_scale*qx_scale*qz*moments[3]*inverse7;
        value[13]=3.0*qx_scale*moments[1]*inverse5-
                  15.0*qx_scale*qy_scale*qy_scale*moments[8]*inverse7;
        value[14]=-15.0*qx_scale*qy_scale*qz*moments[4]*inverse7;
        value[15]=3.0*qx_scale*moments[1]*inverse5-
                  15.0*qx_scale*qz*qz*moments[1]*inverse7;
        value[16]=9.0*qy_scale*moments[2]*inverse5-
                  15.0*qy_scale*qy_scale*qy_scale*moments[9]*inverse7;
        value[17]=3.0*qz*moments[0]*inverse5-
                  15.0*qy_scale*qy_scale*qz*moments[5]*inverse7;
        value[18]=3.0*qy_scale*moments[2]*inverse5-
                  15.0*qy_scale*qz*qz*moments[2]*inverse7;
        value[19]=(9.0*qz*inverse5-15.0*qz*qz*qz*inverse7)*moments[0];

        for(std::size_t component_index=0;component_index<20;++component_index){
            const double term=jacobian*value[component_index];
            compensated_add(term,total[component_index],correction[component_index]);
            magnitude[component_index]+=std::abs(term);
        }
      }
    }
    if(absolute)*absolute=magnitude;
    if(evaluations)*evaluations+=
        radial_panel_count*nodes.size()*nodes.size();
    for(std::size_t component_index=0;component_index<20;++component_index)
        total[component_index]+=correction[component_index];
    return total;
}

inline std::vector<EstimateInfo> polar_cartesian_estimate_all(
    const Prism& prism,const Point& point,const EstimateOptions& options){
    const auto start=std::chrono::steady_clock::now();
    Components previous{},previous_difference{},coarse_difference{};
    Components current{},absolute{};
    std::size_t evaluations=0;
    bool have_previous=false,have_previous_difference=false;
    constexpr std::array<double,5> steps{{
        0.25,0.125,0.0625,0.03125,0.015625}};
    const bool touches_observation_pole=
        (point.lat>=0.0&&close_coordinate(prism.lat2,0.5*PI))||
        (point.lat<0.0&&close_coordinate(prism.lat1,-0.5*PI));
    const double radial_gap=point.r>prism.r2
        ? point.r-prism.r2:prism.r1-point.r;
    const std::size_t step_count=
        touches_observation_pole&&radial_gap/point.r<1.0e-6?5:4;
    for(std::size_t step_index=0;step_index<step_count;++step_index){
        const double step=steps[step_index];
        Components level_absolute{};
        current=polar_cartesian_level(
            prism,point,step,&level_absolute,&evaluations);
        absolute=level_absolute;
        if(have_previous){
            Components difference{};
            for(std::size_t k=0;k<20;++k)
                difference[k]=std::abs(current[k]-previous[k]);
            if(have_previous_difference)
                coarse_difference=previous_difference;
            previous_difference=difference;
            have_previous_difference=true;
        }
        previous=current;
        have_previous=true;
    }

    const double scale=options.gravitational_constant*prism.density;
    const double scale_abs=std::abs(scale);
    std::vector<EstimateInfo> result(20);
    const double milliseconds=std::chrono::duration<double,std::milli>(
        std::chrono::steady_clock::now()-start).count();
    for(std::size_t k=0;k<20;++k){
        auto& item=result[k];
        item.component=static_cast<Component>(k);
        item.value=scale*current[k];
        const double fine_difference=have_previous_difference
            ? previous_difference[k]:0.0;
        const double preceding_difference=coarse_difference[k];
        double normalized_error=fine_difference;
        if(preceding_difference>0.0&&
           fine_difference<0.25*preceding_difference)
            normalized_error=std::min(
                normalized_error,
                2.0*fine_difference*fine_difference/preceding_difference);
        item.error=scale_abs*normalized_error;
        item.absolute=scale_abs*absolute[k];
        item.integral_evaluations=evaluations;
        item.intervals=polar_de_nodes(steps[step_count-1]).size();
        const double effective_rtol=std::max(
            options.integration.rtol,
            500.0*std::numeric_limits<double>::epsilon());
        const double tolerance=std::max({
            options.integration.atol,
            effective_rtol*std::abs(item.value),
            50.0*std::numeric_limits<double>::epsilon()*item.absolute});
        item.converged=std::isfinite(item.value)&&std::isfinite(item.error)&&
                       item.error<=tolerance;
        item.milliseconds=milliseconds;
    }
    return result;
}

} // namespace detail

inline std::vector<EstimateInfo> estimate_all(const Prism& prism,const Point& point,
                                               const EstimateOptions& options={});

inline double estimate_component_independent(Component c,const Prism& prism,
                                              const Point& point,
                                              const EstimateOptions& options,
                                              EstimateInfo* output);

inline double estimate(Component c, const Prism& prism, const Point& point,
                            const EstimateOptions& options={}, EstimateInfo* output=nullptr) {
    return estimate_component_independent(c,prism,point,options,output);
}

inline double estimate(const std::string& name, const Prism& prism,
                            const Point& point, const EstimateOptions& options={},
                            EstimateInfo* info=nullptr) {
    return estimate(component(name),prism,point,options,info);
}

inline double estimate_component_independent(Component component_id,const Prism& prism,
                                              const Point& point,
                                              const EstimateOptions& options,
                                              EstimateInfo* output) {
    const auto start=std::chrono::steady_clock::now();
    detail::validate_geometry(prism,point);
    if(detail::polar_cartesian_applicable(prism,point)){
        auto polar=detail::polar_cartesian_estimate_all(prism,point,options);
        EstimateInfo value=polar[static_cast<std::size_t>(component_id)];
        if(output)*output=value;
        return value.value;
    }
    const detail::BoundaryLocation location=detail::boundary_location(prism,point);
    if(location==detail::BoundaryLocation::EdgeOrCorner)
        throw std::domain_error(
            "observation point lies on a latitude/longitude edge or corner; no boundary convention is defined");
    if(location==detail::BoundaryLocation::LowerRadialFace)
        throw std::domain_error(
            "observation point lies on the lower material face; specify an inner one-sided convention before evaluation");

    if(location==detail::BoundaryLocation::OuterRadialFace){
        double h=std::min(options.radial_limit_relative_step*point.r,
                          0.1*(prism.r2-prism.r1));
        h=std::max(h,256.0*std::numeric_limits<double>::epsilon()*point.r);
        EstimateInfo best;best.component=component_id;
        const std::size_t refinements=std::max<std::size_t>(1,options.radial_limit_max_refinements);
        for(std::size_t refinement=0;refinement<refinements;++refinement){
            Point p1=point,p2=point,p3=point;
            p1.r=prism.r2+h;p2.r=prism.r2+0.5*h;p3.r=prism.r2+0.25*h;
            EstimateInfo i1,i2,i3;
            const double f1=estimate_component_independent(component_id,prism,p1,options,&i1);
            const double f2=estimate_component_independent(component_id,prism,p2,options,&i2);
            const double f3=estimate_component_independent(component_id,prism,p3,options,&i3);
            const double r1=2.0*f2-f1,r2=2.0*f3-f2;
            best=i3;best.value=(4.0*r2-r1)/3.0;
            best.error=(8.0*i3.error+6.0*i2.error+i1.error)/3.0+std::abs(best.value-r2);
            best.absolute=(8.0*i3.absolute+6.0*i2.absolute+i1.absolute)/3.0;
            best.integral_evaluations=i1.integral_evaluations+i2.integral_evaluations+i3.integral_evaluations;
            best.intervals=i1.intervals+i2.intervals+i3.intervals;
            best.milliseconds=i1.milliseconds+i2.milliseconds+i3.milliseconds;
            const double effective_rtol=std::max(options.integration.rtol,
                500.0*std::numeric_limits<double>::epsilon());
            const double tolerance=std::max({options.integration.atol,
                effective_rtol*std::abs(best.value),
                50.0*std::numeric_limits<double>::epsilon()*best.absolute});
            best.converged=i1.converged&&i2.converged&&i3.converged&&best.error<=tolerance;
            if(best.converged)break;
            h*=0.5;
        }
        if(output)*output=best;
        return best.value;
    }

    EstimateInfo result;result.component=component_id;
    const double dr=prism.r2-prism.r1,dlat=prism.lat2-prism.lat1,dlon=prism.lon2-prism.lon1;
    const double physical_scale=0.5*options.gravitational_constant*prism.density*dr*dlat*dlon;
    if(physical_scale==0.0){
        result.milliseconds=std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-start).count();
        if(output)*output=result;
        return 0.0;
    }

    const double lon_center=0.5*(prism.lon1+prism.lon2);
    const bool at_pole=std::abs(std::cos(point.lat))<=64.0*std::numeric_limits<double>::epsilon();
    const double effective_lon=at_pole?0.0:point.lon;
    const std::size_t component_index=static_cast<std::size_t>(component_id);
    const bool symmetry_odd=(component_index==2||component_index==5||component_index==8||
        component_index==11||component_index==14||component_index==16||component_index==18);
    const bool exact_zero=symmetry_odd&&
        std::abs(effective_lon-lon_center)<=64.0*std::numeric_limits<double>::epsilon()*
        std::max(1.0,std::abs(lon_center));
    if(exact_zero){
        result.value=0.0;result.error=0.0;result.converged=true;
        result.milliseconds=std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-start).count();
        if(output)*output=result;
        return 0.0;
    }

    double total=scaled_explicit(component_id,prism,point),correction=0.0,total_error=0.0;
    double total_absolute=std::abs(total);
    bool all_converged=true;
    const double normalized_atol=options.integration.atol/std::max(
        std::abs(physical_scale),std::numeric_limits<double>::min());

    for(int family=0;family<2;++family){
        const EdgeKind kind=family==0?EdgeKind::Latitude:EdgeKind::Longitude;
        const double center=family==0?0.5*(prism.lon1+prism.lon2):0.5*(prism.lat1+prism.lat2);
        const double half=family==0?0.5*dlon:0.5*dlat;
        detail::PreparedPairedKernel prepared;
        {
            const detail::ScopedJetContext context(static_cast<std::size_t>(component_id));
            prepared=prepare_paired_kernel(prism,point,kind);
        }
        auto f=[&](double z){return paired_kernel_prepared(
            component_id,prism,prepared,center+half*z);};

        double fixed=0.0,t0=0.0;
        if(family==0){
            fixed=std::abs(point.lat-prism.lat1)<=std::abs(point.lat-prism.lat2)
                ?prism.lat1:prism.lat2;
            t0=at_pole?center:detail::periodic_candidate(point.lon,prism.lon1,prism.lon2,2.0*PI);
        }else{
            const double d1=std::abs(std::remainder(prism.lon1-effective_lon,2.0*PI));
            const double d2=std::abs(std::remainder(prism.lon2-effective_lon,2.0*PI));
            fixed=d1<=d2?prism.lon1:prism.lon2;
            t0=std::atan2(std::sin(point.lat),std::cos(point.lat)*std::cos(fixed-point.lon));
            t0=detail::clamp(t0,prism.lat1,prism.lat2);
        }
        const double z0=detail::clamp((t0-center)/half,-1.0,1.0);
        double source_lat=t0,source_lon=fixed;
        if(family==0){source_lat=fixed;source_lon=t0;}
        const double gamma=detail::angular_separation(point,source_lat,source_lon);
        double radial_gap=0.0;
        if(point.r<prism.r1)radial_gap=prism.r1-point.r;
        else if(point.r>prism.r2)radial_gap=point.r-prism.r2;
        const double core_width=std::sqrt(gamma*gamma+
            (radial_gap/point.r)*(radial_gap/point.r))/std::max(half,1e-30);
        const bool use_gk21=core_width<1.5&&radial_gap/point.r<0.02;
        const double ztol=256.0*std::numeric_limits<double>::epsilon();
        std::vector<double> split{-1.0,z0,1.0};
        std::sort(split.begin(),split.end());
        split.erase(std::unique(split.begin(),split.end(),[&](double a,double b){
            return std::abs(a-b)<=ztol;}),split.end());
        const std::size_t pieces=split.size()-1;
        for(std::size_t p=0;p<pieces;++p){
            const double a=split[p],b=split[p+1];
            const bool projection_left=std::abs(a-z0)<=ztol;
            const bool projection_right=std::abs(b-z0)<=ztol;
            auto mapped=[&](double s){
                double z=s,jacobian=1.0;
                if(projection_left&&projection_right){
                    z=a+0.5*(b-a)*(1.0-std::cos(s));
                    jacobian=0.5*(b-a)*std::sin(s);
                }else if(projection_left){
                    z=a+(b-a)*(1.0-std::cos(s));
                    jacobian=(b-a)*std::sin(s);
                }else if(projection_right){
                    z=a+(b-a)*std::sin(s);
                    jacobian=(b-a)*std::cos(s);
                }
                return f(z)*jacobian;
            };
            const double integration_a=(projection_left||projection_right)
                ?0.0:a;
            const double integration_b=(projection_left||projection_right)
                ?((projection_left&&projection_right)?PI:0.5*PI):b;
            IntegralOptions piece_options=options.integration;
            piece_options.use_gk21=use_gk21;
            piece_options.atol=normalized_atol/(2.0*static_cast<double>(pieces));
            piece_options.max_evaluations=std::min<std::size_t>(options.integration.max_evaluations,10000);
            piece_options.max_intervals=std::min<std::size_t>(options.integration.max_intervals,500);
            IntegralInfo info;
            const double value=integral(mapped,integration_a,integration_b,
                                        piece_options,&info);
            detail::compensated_add(0.5*value,total,correction);
            total_error+=0.5*info.error;
            total_absolute+=0.5*info.absolute;
            result.integral_evaluations+=info.evaluations;result.intervals+=info.intervals;
            all_converged=all_converged&&info.converged;
        }
    }

    result.value=physical_scale*(total+correction);
    result.error=std::abs(physical_scale)*total_error;
    result.absolute=std::abs(physical_scale)*total_absolute;
    const double effective_rtol=std::max(options.integration.rtol,
        500.0*std::numeric_limits<double>::epsilon());
    const double tolerance=std::max({options.integration.atol,
        effective_rtol*std::abs(result.value),
        50.0*std::numeric_limits<double>::epsilon()*result.absolute});
    result.converged=all_converged&&std::isfinite(result.value)&&
        std::isfinite(result.error)&&result.error<=tolerance;
    result.milliseconds=std::chrono::duration<double,std::milli>(
        std::chrono::steady_clock::now()-start).count();
    if(output)*output=result;
    return result.value;
}

inline std::vector<EstimateInfo> estimate_all(const Prism& prism,const Point& point,
                                               const EstimateOptions& options) {
    detail::validate_geometry(prism,point);
    if(detail::polar_cartesian_applicable(prism,point))
        return detail::polar_cartesian_estimate_all(prism,point,options);
    std::vector<EstimateInfo> result(20);
    for(int k=0;k<20;++k)
        estimate_component_independent(static_cast<Component>(k),prism,point,options,
            &result[static_cast<std::size_t>(k)]);
    return result;
}

// Compatibility interface for the Work1 programs in the original Tesseroid
// package. Input angles are degrees, radii are kilometres, and G=1, exactly as
// in Work1_SetGravityParameters.m. TLI is length-unit invariant, so retaining
// kilometres also retains the original Work output units. Output component
// order is byte-for-byte compatible with Result_TFG.dat.
namespace work_detail {

using OutputPointers = std::array<double*,20>;

inline EstimateOptions options(double atol,double rtol){
    EstimateOptions value;
    value.integration.atol=std::max(0.0,atol);
    value.integration.rtol=std::max(0.0,rtol);
    value.gravitational_constant=1.0;
    return value;
}

inline std::array<double,20> coordinate(
    int model_count,const double* fai1,const double* fai2,
    const double* lamda1,const double* lamda2,const double* r1,
    const double* r2,const double* rho,double latitude,double longitude,
    double radius,const EstimateOptions& opts){
    constexpr double degree=PI/180.0;
    const Point point{degree*latitude,detail::normalize_longitude(degree*longitude),radius};
    std::array<double,20> sum{},correction{};
    for(int j=0;j<model_count;++j){
        const Prism prism{degree*fai1[j],degree*fai2[j],
            degree*lamda1[j],degree*lamda2[j],r1[j],r2[j],rho[j]};
        const auto value=estimate_all(prism,point,opts);
        for(std::size_t k=0;k<sum.size();++k)
            detail::compensated_add(value[k].value,sum[k],correction[k]);
    }
    for(std::size_t k=0;k<sum.size();++k)sum[k]+=correction[k];
    return sum;
}

inline void store(const std::array<double,20>& value,
                  const OutputPointers& output,int i){
    // New canonical order: Vxz/Vyz/Vxyy/Vxzz/Vyzz. The Work file uses the
    // equivalent symmetric labels Vzx/Vzy/Vyyx/Vzzx/Vzzy.
    for(std::size_t k=0;k<value.size();++k)output[k][i]=static_cast<double>(value[k]);
}

using MagneticOutputPointers = std::array<double*,10>;

inline std::array<double,10> magnetic_coordinate(
    int model_count,const double* fai1,const double* fai2,
    const double* lamda1,const double* lamda2,const double* r1,
    const double* r2,const double* mx,const double* my,const double* mz,
    double latitude,double longitude,double radius,const EstimateOptions& opts){
    constexpr double degree=PI/180.0;
    const Point point{degree*latitude,detail::normalize_longitude(degree*longitude),radius};
    const double data_lat=degree*latitude;
    const double cos_data_lat=std::cos(data_lat);
    const double sin_data_lat=std::sin(data_lat);
    std::array<double,10> sum{},correction{};
    for(int j=0;j<model_count;++j){
        const Prism prism{degree*fai1[j],degree*fai2[j],
            degree*lamda1[j],degree*lamda2[j],r1[j],r2[j],1.0};
        const auto g=estimate_all(prism,point,opts);

        // Rotate the model-centred north/east/up magnetization into the
        // observation-centred Cartesian frame. This is the same rotation as
        // the original Tesseroid_MagneticEstimate implementation.
        const double beta=degree*(0.5*(lamda1[j]+lamda2[j])-longitude);
        const double source_lat=degree*0.5*(fai1[j]+fai2[j]);
        const double cb=std::cos(beta),sb=std::sin(beta);
        const double cs=std::cos(source_lat),ss=std::sin(source_lat);
        const double ml_x=(cb*sin_data_lat*ss+cos_data_lat*cs)*mx[j]+
            (sb*sin_data_lat)*my[j]-(cb*sin_data_lat*cs-cos_data_lat*ss)*mz[j];
        const double ml_y=(-sb*ss)*mx[j]+cb*my[j]+sb*cs*mz[j];
        const double ml_z=-(cb*cos_data_lat*ss-sin_data_lat*cs)*mx[j]-
            (sb*cos_data_lat)*my[j]+(cb*cos_data_lat*cs+sin_data_lat*ss)*mz[j];
        auto value=[&](std::size_t k){return g[k].value;};

        const std::array<double,10> contribution{
            ml_x*value(1)+ml_y*value(2)+ml_z*value(3),
            ml_x*value(4)+ml_y*value(5)+ml_z*value(6),
            ml_x*value(5)+ml_y*value(7)+ml_z*value(8),
            ml_x*value(6)+ml_y*value(8)+ml_z*value(9),
            ml_x*value(10)+ml_y*value(11)+ml_z*value(12),
            ml_x*value(11)+ml_y*value(13)+ml_z*value(14),
            ml_x*value(13)+ml_y*value(16)+ml_z*value(17),
            ml_x*value(12)+ml_y*value(14)+ml_z*value(15),
            ml_x*value(14)+ml_y*value(17)+ml_z*value(18),
            ml_x*value(15)+ml_y*value(18)+ml_z*value(19)};
        for(std::size_t k=0;k<sum.size();++k)
            detail::compensated_add(contribution[k],sum[k],correction[k]);
    }
    for(std::size_t k=0;k<sum.size();++k)sum[k]+=correction[k];
    return sum;
}

inline void store(const std::array<double,10>& value,
                  const MagneticOutputPointers& output,int i){
    for(std::size_t k=0;k<value.size();++k)output[k][i]=static_cast<double>(value[k]);
}

} // namespace work_detail

inline void Tesseroid_GravityEstimate(
    int NumDensityModel,const double* Fai1,const double* Fai2,
    const double* Lamda1,const double* Lamda2,const double* R1,
    const double* R2,const double* rho,int NumCoordinates,
    const double* Latitude,const double* Longitude,const double* Radius,
    double AbsTol,double RelTol,double* V,double* Vx,double* Vy,double* Vz,
    double* Vxx,double* Vxy,double* Vyy,double* Vzx,double* Vzy,double* Vzz,
    double* Vxxx,double* Vxxy,double* Vxxz,double* Vxyz,double* Vyyx,
    double* Vyyy,double* Vyyz,double* Vzzx,double* Vzzy,double* Vzzz){
    const auto opts=work_detail::options(AbsTol,RelTol);
    const work_detail::OutputPointers output{V,Vx,Vy,Vz,Vxx,Vxy,Vzx,Vyy,Vzy,Vzz,
        Vxxx,Vxxy,Vxxz,Vyyx,Vxyz,Vzzx,Vyyy,Vyyz,Vzzy,Vzzz};
    for(int i=0;i<NumCoordinates;++i){
        const auto value=work_detail::coordinate(NumDensityModel,Fai1,Fai2,Lamda1,
            Lamda2,R1,R2,rho,Latitude[i],Longitude[i],Radius[i],opts);
        work_detail::store(value,output,i);
    }
}

inline void Tesseroid_GravityEstimateOpenMP(
    int NumDensityModel,const double* Fai1,const double* Fai2,
    const double* Lamda1,const double* Lamda2,const double* R1,
    const double* R2,const double* rho,int NumCoordinates,
    const double* Latitude,const double* Longitude,const double* Radius,
    double AbsTol,double RelTol,double* V,double* Vx,double* Vy,double* Vz,
    double* Vxx,double* Vxy,double* Vyy,double* Vzx,double* Vzy,double* Vzz,
    double* Vxxx,double* Vxxy,double* Vxxz,double* Vxyz,double* Vyyx,
    double* Vyyy,double* Vyyz,double* Vzzx,double* Vzzy,double* Vzzz){
    const auto opts=work_detail::options(AbsTol,RelTol);
    const work_detail::OutputPointers output{V,Vx,Vy,Vz,Vxx,Vxy,Vzx,Vyy,Vzy,Vzz,
        Vxxx,Vxxy,Vxxz,Vyyx,Vxyz,Vzzx,Vyyy,Vyyz,Vzzy,Vzzz};
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i=0;i<NumCoordinates;++i){
        const auto value=work_detail::coordinate(NumDensityModel,Fai1,Fai2,Lamda1,
            Lamda2,R1,R2,rho,Latitude[i],Longitude[i],Radius[i],opts);
        work_detail::store(value,output,i);
    }
}

inline void Tesseroid_MagneticEstimate(
    int NumMagnetizationModel,const double* Fai1,const double* Fai2,
    const double* Lamda1,const double* Lamda2,const double* R1,
    const double* R2,const double* Mx,const double* My,const double* Mz,
    int NumCoordinates,const double* Latitude,const double* Longitude,
    const double* Radius,double AbsTol,double RelTol,double* V,double* Vx,
    double* Vy,double* Vz,double* Vxx,double* Vxy,double* Vyy,double* Vzx,
    double* Vzy,double* Vzz){
    const auto opts=work_detail::options(AbsTol,RelTol);
    const work_detail::MagneticOutputPointers output{V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz};
    for(int i=0;i<NumCoordinates;++i){
        const auto value=work_detail::magnetic_coordinate(NumMagnetizationModel,
            Fai1,Fai2,Lamda1,Lamda2,R1,R2,Mx,My,Mz,Latitude[i],Longitude[i],
            Radius[i],opts);
        work_detail::store(value,output,i);
    }
}

inline void Tesseroid_MagneticEstimateOpenMP(
    int NumMagnetizationModel,const double* Fai1,const double* Fai2,
    const double* Lamda1,const double* Lamda2,const double* R1,
    const double* R2,const double* Mx,const double* My,const double* Mz,
    int NumCoordinates,const double* Latitude,const double* Longitude,
    const double* Radius,double AbsTol,double RelTol,double* V,double* Vx,
    double* Vy,double* Vz,double* Vxx,double* Vxy,double* Vyy,double* Vzx,
    double* Vzy,double* Vzz){
    const auto opts=work_detail::options(AbsTol,RelTol);
    const work_detail::MagneticOutputPointers output{V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz};
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i=0;i<NumCoordinates;++i){
        const auto value=work_detail::magnetic_coordinate(NumMagnetizationModel,
            Fai1,Fai2,Lamda1,Lamda2,R1,R2,Mx,My,Mz,Latitude[i],Longitude[i],
            Radius[i],opts);
        work_detail::store(value,output,i);
    }
}

#ifdef TLI_ENABLE_MPI
inline void Tesseroid_GravityEstimateMPI(
    int NumDensityModel,const double* Fai1,const double* Fai2,
    const double* Lamda1,const double* Lamda2,const double* R1,
    const double* R2,const double* rho,int NumCoordinates,
    const double* Latitude,const double* Longitude,const double* Radius,
    double AbsTol,double RelTol,double* V,double* Vx,double* Vy,double* Vz,
    double* Vxx,double* Vxy,double* Vyy,double* Vzx,double* Vzy,double* Vzz,
    double* Vxxx,double* Vxxy,double* Vxxz,double* Vxyz,double* Vyyx,
    double* Vyyy,double* Vyyz,double* Vzzx,double* Vzzy,double* Vzzz){
    int rank=0,size=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    std::vector<int> counts(size),displacements(size);
    for(int p=0;p<size;++p){
        counts[p]=NumCoordinates/size+(p<NumCoordinates%size?1:0);
        displacements[p]=p==0?0:displacements[p-1]+counts[p-1];
    }
    const int first=displacements[rank],local_count=counts[rank];
    std::array<std::vector<double>,20> local;
    for(auto& component:local)component.resize(static_cast<std::size_t>(local_count));
    std::array<double*,20> local_output{};
    for(std::size_t k=0;k<local.size();++k)local_output[k]=local[k].data();
    const auto opts=work_detail::options(AbsTol,RelTol);
    for(int j=0;j<local_count;++j){
        const int i=first+j;
        const auto value=work_detail::coordinate(NumDensityModel,Fai1,Fai2,Lamda1,
            Lamda2,R1,R2,rho,Latitude[i],Longitude[i],Radius[i],opts);
        work_detail::store(value,local_output,j);
    }
    const work_detail::OutputPointers output{V,Vx,Vy,Vz,Vxx,Vxy,Vzx,Vyy,Vzy,Vzz,
        Vxxx,Vxxy,Vxxz,Vyyx,Vxyz,Vzzx,Vyyy,Vyyz,Vzzy,Vzzz};
    for(std::size_t k=0;k<local.size();++k)
        MPI_Gatherv(local[k].data(),local_count,MPI_DOUBLE,output[k],counts.data(),
            displacements.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

inline void Tesseroid_MagneticEstimateMPI(
    int NumMagnetizationModel,const double* Fai1,const double* Fai2,
    const double* Lamda1,const double* Lamda2,const double* R1,
    const double* R2,const double* Mx,const double* My,const double* Mz,
    int NumCoordinates,const double* Latitude,const double* Longitude,
    const double* Radius,double AbsTol,double RelTol,double* V,double* Vx,
    double* Vy,double* Vz,double* Vxx,double* Vxy,double* Vyy,double* Vzx,
    double* Vzy,double* Vzz){
    int rank=0,size=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    std::vector<int> counts(size),displacements(size);
    for(int p=0;p<size;++p){
        counts[p]=NumCoordinates/size+(p<NumCoordinates%size?1:0);
        displacements[p]=p==0?0:displacements[p-1]+counts[p-1];
    }
    const int first=displacements[rank],local_count=counts[rank];
    std::array<std::vector<double>,10> local;
    for(auto& component:local)component.resize(static_cast<std::size_t>(local_count));
    std::array<double*,10> local_output{};
    for(std::size_t k=0;k<local.size();++k)local_output[k]=local[k].data();
    const auto opts=work_detail::options(AbsTol,RelTol);
    for(int j=0;j<local_count;++j){
        const int i=first+j;
        const auto value=work_detail::magnetic_coordinate(NumMagnetizationModel,
            Fai1,Fai2,Lamda1,Lamda2,R1,R2,Mx,My,Mz,Latitude[i],Longitude[i],
            Radius[i],opts);
        work_detail::store(value,local_output,j);
    }
    const work_detail::MagneticOutputPointers output{V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz};
    for(std::size_t k=0;k<local.size();++k)
        MPI_Gatherv(local[k].data(),local_count,MPI_DOUBLE,output[k],counts.data(),
            displacements.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
}
#endif

} // namespace Tesseroid
