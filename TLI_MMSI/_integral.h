#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <vector>

namespace Tesseroid {

struct IntegralInfo {
    double error = 0.0;
    // Approximation to the integral of |f| used by the global roundoff floor.
    double absolute = 0.0;
    std::size_t evaluations = 0;
    std::size_t intervals = 0;
    bool converged = true;
};

struct IntegralOptions {
    double atol = 1.0e-12;
    double rtol = 1.0e-10;
    std::size_t max_evaluations = 200000;
    std::size_t max_intervals = 10000;
    // GK21 is advantageous for nearby or structured kernels; callers may
    // select it per integral while retaining GK15 for inexpensive far fields.
    bool use_gk21 = false;
};

namespace detail {

inline void compensated_add(double value,double& sum,double& correction){
    const double updated=sum+value;
    correction+=(std::abs(sum)>=std::abs(value))
        ? (sum-updated)+value
        : (value-updated)+sum;
    sum=updated;
}

struct GKInterval {
    double a = 0.0;
    double b = 0.0;
    double value = 0.0;
    double error = 0.0;
    double absolute = 0.0;
};

struct LargerError {
    bool operator()(const GKInterval& lhs, const GKInterval& rhs) const {
        return lhs.error < rhs.error;
    }
};

template <class F>
inline GKInterval gk15(const F& f, double a, double b,
                       std::size_t& evaluations) {
    static constexpr std::array<double, 8> xgk = {
        0.991455371120812639206854697526329,
        0.949107912342758524526189684047851,
        0.864864423359769072789712788640926,
        0.741531185599394439863864773280788,
        0.586087235467691130294144838258730,
        0.405845151377397166906606412076961,
        0.207784955007898467600689403773245,
        0.0};
    static constexpr std::array<double, 4> wg = {
        0.129484966168869693270611432679082,
        0.279705391489276667901467771423780,
        0.381830050505118944950369775488975,
        0.417959183673469387755102040816327};
    static constexpr std::array<double, 8> wgk = {
        0.022935322010529224963732008058970,
        0.063092092629978553290700663189204,
        0.104790010322250183839876322541518,
        0.140653259715525918745189590510238,
        0.169004726639267902826583426598550,
        0.190350578064785409913256402421014,
        0.204432940075298892414161999234649,
        0.209482141084727828012999174891714};

    const double center = 0.5 * (a + b);
    const double half = 0.5 * (b - a);
    const double fc = f(center);
    evaluations += 1;

    double resg = wg[3] * fc;
    double resk = wgk[7] * fc;
    double resabs = wgk[7] * std::abs(fc);
    std::array<double, 7> fsum{}, fleft{}, fright{};

    for (std::size_t j = 0; j < 7; ++j) {
        const double abscissa = half * xgk[j];
        const double f1 = f(center - abscissa);
        const double f2 = f(center + abscissa);
        evaluations += 2;
        fleft[j] = f1;
        fright[j] = f2;
        fsum[j] = f1 + f2;
        resk = std::fma(wgk[j],fsum[j],resk);
        resabs = std::fma(wgk[j],std::abs(f1)+std::abs(f2),resabs);
    }
    resg += wg[0] * fsum[1] + wg[1] * fsum[3] + wg[2] * fsum[5];

    const double mean = 0.5 * resk;
    double resasc = wgk[7] * std::abs(fc - mean);
    for (std::size_t j = 0; j < 7; ++j) {
        resasc += wgk[j] * (std::abs(fleft[j] - mean) + std::abs(fright[j] - mean));
    }

    const double value = resk * half;
    resabs *= std::abs(half);
    resasc *= std::abs(half);
    double error = std::abs((resk - resg) * half);
    if (resasc > 0.0 && error > 0.0) {
        const double scale = std::pow(200.0 * error / resasc, 1.5);
        error = resasc * std::min(1.0, scale);
    }
    const double roundoff = 50.0 * std::numeric_limits<double>::epsilon() * resabs;
    error = std::max(error, roundoff);
    return GKInterval{a, b, value, error, resabs};
}

template <class F>
inline GKInterval gk21(const F& f,double a,double b,
                       std::size_t& evaluations){
    static constexpr std::array<double,11> xgk={
        0.995657163025808080735527280689003,0.973906528517171720077964012084452,
        0.930157491355708226001207180059508,0.865063366688984510732096688423493,
        0.780817726586416897063717578345042,0.679409568299024406234327365114874,
        0.562757134668604683339000099272694,0.433395394129247190799265943165784,
        0.294392862701460198131126603103866,0.148874338981631210884826001129720,
        0.0};
    static constexpr std::array<double,5> wg={
        0.066671344308688137593568809893332,0.149451349150580593145776339657697,
        0.219086362515982043995534934228163,0.269266719309996355091226921569469,
        0.295524224714752870173892994651338};
    static constexpr std::array<double,11> wgk={
        0.011694638867371874278064396062192,0.032558162307964727478818972459390,
        0.054755896574351996031381300244580,0.075039674810919952767043140916190,
        0.093125454583697605535065465083366,0.109387158802297641899210590325805,
        0.123491976262065851077958109831074,0.134709217311473325928054001771707,
        0.142775938577060080797094273138717,0.147739104901338491374841515972068,
        0.149445554002916905664936468389821};
    const double center=0.5*(a+b),half=0.5*(b-a);
    const double fc=f(center);++evaluations;
    double resg=0.0,resk=wgk[10]*fc,resabs=wgk[10]*std::abs(fc);
    std::array<double,10> fl{},fr{},fs{};
    for(std::size_t j=0;j<10;++j){
        const double dx=half*xgk[j];
        fl[j]=f(center-dx);fr[j]=f(center+dx);evaluations+=2;
        fs[j]=fl[j]+fr[j];
        resk=std::fma(wgk[j],fs[j],resk);
        resabs=std::fma(wgk[j],std::abs(fl[j])+std::abs(fr[j]),resabs);
        if((j&1U)!=0U)resg=std::fma(wg[(j-1)/2],fs[j],resg);
    }
    const double mean=0.5*resk;
    double resasc=wgk[10]*std::abs(fc-mean);
    for(std::size_t j=0;j<10;++j)
        resasc+=wgk[j]*(std::abs(fl[j]-mean)+std::abs(fr[j]-mean));
    const double value=resk*half;
    resabs*=std::abs(half);resasc*=std::abs(half);
    double error=std::abs((resk-resg)*half);
    if(resasc>0.0&&error>0.0){
        const double scale=std::pow(200.0*error/resasc,1.5);
        error=resasc*std::min(1.0,scale);
    }
    error=std::max(error,50.0*std::numeric_limits<double>::epsilon()*resabs);
    return GKInterval{a,b,value,error,resabs};
}

template<class F>
inline GKInterval scalar_rule(const F& f,double a,double b,
                              std::size_t& evaluations,bool use_gk21){
    return use_gk21?gk21(f,a,b,evaluations):gk15(f,a,b,evaluations);
}

inline std::size_t scalar_rule_cost(bool use_gk21){
    return use_gk21?21U:15U;
}

} // namespace detail

// Adaptive Gauss--Kronrod 7/15.  The accepted global tolerance is
// max(atol, rtol*abs(I), 50*epsilon*integral(abs(f))), matching the
// manuscript's absolute, relative, and attainable-roundoff stopping rule.
template <class F>
inline double integral(const F& f, double a, double b,
                            const IntegralOptions& options = {},
                            IntegralInfo* info = nullptr) {
    IntegralInfo local;
    if (a == b) {
        if (info) *info = local;
        return 0.0;
    }
    if (b < a) {
        const double value = -integral(f, b, a, options, &local);
        if (info) *info = local;
        return value;
    }

    std::priority_queue<detail::GKInterval, std::vector<detail::GKInterval>,
                        detail::LargerError> queue;
    auto first=detail::scalar_rule(f,a,b,local.evaluations,options.use_gk21);
    queue.push(first);
    double total = first.value;
    double correction = 0.0;
    double total_error = first.error;
    double total_absolute = first.absolute;
    local.intervals = 1;

    const double rtol = std::max(options.rtol,
        500.0 * std::numeric_limits<double>::epsilon());
    auto tolerance = [&]() {
        const double roundoff_floor=50.0*std::numeric_limits<double>::epsilon()*
            total_absolute;
        return std::max({std::max(0.0,options.atol),
            rtol*std::abs(total+correction),roundoff_floor});
    };

    while (total_error > tolerance()) {
        if (queue.empty() || local.evaluations+2*detail::scalar_rule_cost(options.use_gk21)>options.max_evaluations ||
            local.intervals + 1 > options.max_intervals) {
            local.converged = false;
            break;
        }
        const auto worst = queue.top();
        queue.pop();
        const double mid = 0.5 * (worst.a + worst.b);
        if (mid == worst.a || mid == worst.b) {
            local.converged = false;
            break;
        }
        auto left=detail::scalar_rule(f,worst.a,mid,local.evaluations,options.use_gk21);
        auto right=detail::scalar_rule(f,mid,worst.b,local.evaluations,options.use_gk21);
        const double delta=(left.value-worst.value)+right.value;
        detail::compensated_add(delta,total,correction);
        total_error += left.error + right.error - worst.error;
        total_error = std::max(0.0, total_error);
        total_absolute += left.absolute + right.absolute - worst.absolute;
        total_absolute = std::max(0.0,total_absolute);
        queue.push(left);
        queue.push(right);
        ++local.intervals;
    }

    local.error = total_error;
    local.absolute = total_absolute;
    if (info) *info = local;
    return total+correction;
}

} // namespace Tesseroid
