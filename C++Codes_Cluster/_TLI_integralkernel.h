#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>


namespace Tesseroid {

constexpr double PI = 3.141592653589793238462643383279502884;

enum class Component {
    V, Vx, Vy, Vz,
    Vxx, Vxy, Vxz, Vyy, Vyz, Vzz,
    Vxxx, Vxxy, Vxxz, Vxyy, Vxyz, Vxzz,
    Vyyy, Vyyz, Vyzz, Vzzz
};

inline const char* component_name(Component c) {
    static constexpr const char* names[] = {
        "V", "Vx", "Vy", "Vz", "Vxx", "Vxy", "Vxz", "Vyy", "Vyz", "Vzz",
        "Vxxx", "Vxxy", "Vxxz", "Vxyy", "Vxyz", "Vxzz", "Vyyy", "Vyyz", "Vyzz", "Vzzz"};
    return names[static_cast<std::size_t>(c)];
}

inline Component component(const std::string& name) {
    for (int i = 0; i < 20; ++i) {
        const auto c = static_cast<Component>(i);
        if (name == component_name(c)) return c;
    }
    if (name == "Vzx") return Component::Vxz;
    if (name == "Vzy") return Component::Vyz;
    if (name == "Vyyx") return Component::Vxyy;
    if (name == "Vzzx") return Component::Vxzz;
    if (name == "Vzzy") return Component::Vyzz;
    throw std::invalid_argument("unknown gravity component: " + name);
}

inline int component_order(Component c) {
    const int i = static_cast<int>(c);
    return i == 0 ? 0 : (i <= 3 ? 1 : (i <= 9 ? 2 : 3));
}

struct Prism {
    double lat1 = 0.0;
    double lat2 = 0.0;
    double lon1 = 0.0;
    double lon2 = 0.0;
    double r1 = 0.0;
    double r2 = 0.0;
    double density = 0.0;
};

struct Point {
    double lat = 0.0;
    double lon = 0.0;
    double r = 0.0;
};

enum class EdgeKind { Latitude, Longitude };

struct Edge {
    EdgeKind kind = EdgeKind::Latitude;
    double fixed = 0.0;
    double a = 0.0;
    double b = 0.0;
    int sign = 1;
};

namespace detail {

// Compact graded-lexicographic index for x+y+z<=3.  The previous 4x4x4
// storage carried 44 coefficients that are mathematically outside a
// third-order jet; eliminating them reduces every Jet copy, scale, and zero
// initialization without changing any retained coefficient.
inline constexpr std::array<int,64> jet_index{{
     0, 3, 9,19,  2, 8,18, 0,  7,17, 0, 0, 16, 0, 0, 0,
     1, 6,15, 0,  5,14, 0, 0, 13, 0, 0, 0,  0, 0, 0, 0,
     4,12, 0, 0, 11, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    10, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0}};
inline constexpr int ji(int x,int y,int z) { return jet_index[(x*4+y)*4+z]; }

struct JetExponent { std::uint8_t x,y,z; };
inline constexpr std::array<JetExponent,20> jet_exponent{{
    {0,0,0},{1,0,0},{0,1,0},{0,0,1},
    {2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2},
    {3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},{1,0,2},
    {0,3,0},{0,2,1},{0,1,2},{0,0,3}}};
struct JetProduct { std::uint8_t left,right,out; };
inline constexpr std::array<JetProduct,84> make_jet_products(){
    std::array<JetProduct,84> result{};std::size_t n=0;
    for(std::size_t i=0;i<jet_exponent.size();++i)
    for(std::size_t j=0;j<jet_exponent.size();++j){
        const int x=jet_exponent[i].x+jet_exponent[j].x;
        const int y=jet_exponent[i].y+jet_exponent[j].y;
        const int z=jet_exponent[i].z+jet_exponent[j].z;
        if(x+y+z<=3)result[n++]={static_cast<std::uint8_t>(i),
            static_cast<std::uint8_t>(j),static_cast<std::uint8_t>(ji(x,y,z))};
    }
    return result;
}
// Keep table construction compatible with older cluster libstdc++ versions;
// it is initialized once and remains immutable after startup.
inline const auto jet_products=make_jet_products();

struct JetContext {
    std::array<std::uint8_t,20> indices{};
    std::array<JetProduct,84> products{};
    std::size_t index_count=0,product_count=0;
};
inline const std::array<JetContext,21>& jet_context_table(){
    static const std::array<JetContext,21> table=[] {
        std::array<JetContext,21> result{};
        for(std::size_t target=0;target<result.size();++target){
            const bool full=target==20;
            const auto alpha=full?JetExponent{3,3,3}:jet_exponent[target];
            auto& context=result[target];
            for(std::size_t i=0;i<jet_exponent.size();++i){
                const auto beta=jet_exponent[i];
                if(full||(beta.x<=alpha.x&&beta.y<=alpha.y&&beta.z<=alpha.z))
                    context.indices[context.index_count++]=static_cast<std::uint8_t>(i);
            }
            for(const auto& product:jet_products){
                const auto beta=jet_exponent[product.out];
                if(full||(beta.x<=alpha.x&&beta.y<=alpha.y&&beta.z<=alpha.z))
                    context.products[context.product_count++]=product;
            }
        }
        return result;
    }();
    return table;
}
inline thread_local const JetContext* active_jet_context=&jet_context_table()[20];
struct ScopedJetContext {
    const JetContext* previous;
    explicit ScopedJetContext(std::size_t component_index):previous(active_jet_context){
        active_jet_context=&jet_context_table()[component_index];
    }
    ~ScopedJetContext(){active_jet_context=previous;}
};

// Normalized multivariate Taylor coefficients through total degree three.
struct Jet {
    std::array<double, 20> a{};

    Jet() = default;
    Jet(double value) { a[0] = value; }

    double value() const { return a[0]; }
    double& at(int x, int y, int z) { return a[ji(x, y, z)]; }
    double at(int x, int y, int z) const { return a[ji(x, y, z)]; }

    static Jet variable(double value, int axis) {
        Jet result(value);
        if (axis == 0) result.at(1, 0, 0) = 1.0;
        if (axis == 1) result.at(0, 1, 0) = 1.0;
        if (axis == 2) result.at(0, 0, 1) = 1.0;
        return result;
    }
};

inline Jet operator+(const Jet& x, const Jet& y) {
    Jet z;
    for(std::size_t n=0;n<active_jet_context->index_count;++n){
        const auto i=active_jet_context->indices[n];z.a[i]=x.a[i]+y.a[i];
    }
    return z;
}
inline Jet operator-(const Jet& x, const Jet& y) {
    Jet z;
    for(std::size_t n=0;n<active_jet_context->index_count;++n){
        const auto i=active_jet_context->indices[n];z.a[i]=x.a[i]-y.a[i];
    }
    return z;
}
inline Jet operator-(const Jet& x) {
    Jet z;
    for(std::size_t n=0;n<active_jet_context->index_count;++n){
        const auto i=active_jet_context->indices[n];z.a[i]=-x.a[i];
    }
    return z;
}
inline Jet operator*(const Jet& x, const Jet& y) {
    Jet z;
    for(std::size_t n=0;n<active_jet_context->product_count;++n){
        const auto& term=active_jet_context->products[n];
        z.a[term.out]+=x.a[term.left]*y.a[term.right];
    }
    return z;
}
inline Jet operator*(const Jet& x, double y) {
    Jet z;
    for(std::size_t n=0;n<active_jet_context->index_count;++n){
        const auto i=active_jet_context->indices[n];z.a[i]=x.a[i]*y;
    }
    return z;
}
inline Jet operator*(double x, const Jet& y) { return y * x; }
inline Jet operator/(const Jet& x, double y) { return x * (1.0 / y); }
inline Jet operator+(const Jet& x, double y) { Jet z = x; z.a[0] += y; return z; }
inline Jet operator+(double x, const Jet& y) { return y + x; }
inline Jet operator-(const Jet& x, double y) { Jet z = x; z.a[0] -= y; return z; }
inline Jet operator-(double x, const Jet& y) { return Jet(x) - y; }

inline Jet compose3(const Jet& x, double f0, double f1,
                    double f2, double f3) {
    Jet h = x - x.value();
    const Jet h2 = h * h;
    return Jet(f0) + f1 * h + f2 * h2 + f3 * (h2 * h);
}

inline Jet inverse(const Jet& x) {
    const double v = x.value();
    return compose3(x, 1.0 / v, -1.0 / (v * v),
                    1.0 / (v * v * v), -1.0 / (v * v * v * v));
}
inline Jet operator/(const Jet& x, const Jet& y) { return x * inverse(y); }
inline Jet sqrt(const Jet& x) {
    const double v = x.value();
    const double s = std::sqrt(v);
    return compose3(x, s, 1.0 / (2.0 * s),
                    -1.0 / (8.0 * v * s), 1.0 / (16.0 * v * v * s));
}
inline Jet log(const Jet& x) {
    const double v = x.value();
    return compose3(x, std::log(v), 1.0 / v,
                    -1.0 / (2.0 * v * v), 1.0 / (3.0 * v * v * v));
}

// Auxiliary normalized Taylor jet in one source parameter.  Each coefficient
// is itself an observation-coordinate Jet, so differentiating a local divided
// difference cannot reintroduce endpoint cancellation (Test.tex, Eq. 5.5.1).
template<std::size_t N> struct SourceJet {
    std::array<Jet,N+1> a{};
    SourceJet() = default;
    SourceJet(const Jet& value) { a[0]=value; }
    static SourceJet variable(const Jet& value) {
        SourceJet z(value);z.a[1]=Jet(1.0);return z;
    }
};

template<std::size_t N> inline SourceJet<N> operator+(const SourceJet<N>&x,const SourceJet<N>&y){
    SourceJet<N>z;for(std::size_t i=0;i<=N;++i)z.a[i]=x.a[i]+y.a[i];return z;
}
template<std::size_t N> inline SourceJet<N> operator-(const SourceJet<N>&x,const SourceJet<N>&y){
    SourceJet<N>z;for(std::size_t i=0;i<=N;++i)z.a[i]=x.a[i]-y.a[i];return z;
}
template<std::size_t N> inline SourceJet<N> operator-(const SourceJet<N>&x){
    SourceJet<N>z;for(std::size_t i=0;i<=N;++i)z.a[i]=-x.a[i];return z;
}
template<std::size_t N> inline SourceJet<N> operator+(const SourceJet<N>&x,const Jet&y){SourceJet<N>z=x;z.a[0]=z.a[0]+y;return z;}
template<std::size_t N> inline SourceJet<N> operator+(const Jet&x,const SourceJet<N>&y){return y+x;}
template<std::size_t N> inline SourceJet<N> operator-(const SourceJet<N>&x,const Jet&y){SourceJet<N>z=x;z.a[0]=z.a[0]-y;return z;}
template<std::size_t N> inline SourceJet<N> operator-(const Jet&x,const SourceJet<N>&y){return SourceJet<N>{x}-y;}
template<std::size_t N> inline SourceJet<N> operator*(const SourceJet<N>&x,const SourceJet<N>&y){
    SourceJet<N>z;for(std::size_t n=0;n<=N;++n)for(std::size_t k=0;k<=n;++k)z.a[n]=z.a[n]+x.a[k]*y.a[n-k];return z;
}
template<std::size_t N> inline SourceJet<N> operator*(const SourceJet<N>&x,const Jet&y){
    SourceJet<N>z;for(std::size_t i=0;i<=N;++i)z.a[i]=x.a[i]*y;return z;
}
template<std::size_t N> inline SourceJet<N> operator*(const Jet&x,const SourceJet<N>&y){return y*x;}
template<std::size_t N> inline SourceJet<N> inverse(const SourceJet<N>&x){
    SourceJet<N>z;z.a[0]=inverse(x.a[0]);
    for(std::size_t n=1;n<=N;++n){Jet s;for(std::size_t k=1;k<=n;++k)s=s+x.a[k]*z.a[n-k];z.a[n]=-s*z.a[0];}
    return z;
}
template<std::size_t N> inline SourceJet<N> operator/(const SourceJet<N>&x,const SourceJet<N>&y){return x*inverse(y);}
template<std::size_t N> inline SourceJet<N> sqrt(const SourceJet<N>&x){
    SourceJet<N>z;z.a[0]=sqrt(x.a[0]);
    for(std::size_t n=1;n<=N;++n){Jet s;for(std::size_t k=1;k<n;++k)s=s+z.a[k]*z.a[n-k];z.a[n]=(x.a[n]-s)/(2.0*z.a[0]);}
    return z;
}
template<std::size_t N> inline SourceJet<N> log(const SourceJet<N>&x){
    SourceJet<N>xp;for(std::size_t n=0;n<N;++n)xp.a[n]=x.a[n+1]*static_cast<double>(n+1);
    const SourceJet<N>q=xp/x;SourceJet<N>z;z.a[0]=log(x.a[0]);
    for(std::size_t n=0;n<N;++n)z.a[n+1]=q.a[n]/static_cast<double>(n+1);return z;
}

template<std::size_t N> struct BiSourceJet {
    std::array<std::array<Jet,N+1>,N+1> a{};
    BiSourceJet()=default;BiSourceJet(const Jet&v){a[0][0]=v;}
    static BiSourceJet variable_a(const Jet&v){BiSourceJet z(v);z.a[1][0]=Jet(1);return z;}
    static BiSourceJet variable_c(const Jet&v){BiSourceJet z(v);z.a[0][1]=Jet(1);return z;}
};
template<std::size_t N> inline BiSourceJet<N> operator+(const BiSourceJet<N>&x,const BiSourceJet<N>&y){BiSourceJet<N>z;for(std::size_t i=0;i<=N;++i)for(std::size_t j=0;j<=N;++j)z.a[i][j]=x.a[i][j]+y.a[i][j];return z;}
template<std::size_t N> inline BiSourceJet<N> operator-(const BiSourceJet<N>&x,const BiSourceJet<N>&y){BiSourceJet<N>z;for(std::size_t i=0;i<=N;++i)for(std::size_t j=0;j<=N;++j)z.a[i][j]=x.a[i][j]-y.a[i][j];return z;}
template<std::size_t N> inline BiSourceJet<N> operator-(const BiSourceJet<N>&x){BiSourceJet<N>z;for(std::size_t i=0;i<=N;++i)for(std::size_t j=0;j<=N;++j)z.a[i][j]=-x.a[i][j];return z;}
template<std::size_t N> inline BiSourceJet<N> operator*(const BiSourceJet<N>&x,const BiSourceJet<N>&y){
    BiSourceJet<N> z;
    for(std::size_t i=0;i<=N;++i)
        for(std::size_t j=0;j<=N;++j)
            for(std::size_t p=0;p<=i;++p)
                for(std::size_t q=0;q<=j;++q)
                    z.a[i][j]=z.a[i][j]+x.a[p][q]*y.a[i-p][j-q];
    return z;
}
template<std::size_t N> inline BiSourceJet<N> operator*(const BiSourceJet<N>&x,const Jet&y){BiSourceJet<N>z;for(std::size_t i=0;i<=N;++i)for(std::size_t j=0;j<=N;++j)z.a[i][j]=x.a[i][j]*y;return z;}
template<std::size_t N> inline BiSourceJet<N> operator*(const Jet&x,const BiSourceJet<N>&y){return y*x;}
template<std::size_t N> inline BiSourceJet<N> bi_inverse(const BiSourceJet<N>&x){
    const Jet x0=x.a[0][0],ix0=inverse(x0);BiSourceJet<N>h=x-BiSourceJet<N>{x0};h=h*ix0;
    BiSourceJet<N>sum{Jet(1)},power{Jet(1)};double sign=-1;
    for(std::size_t k=1;k<=2*N;++k){power=power*h;sum=sum+Jet(sign)*power;sign=-sign;}return ix0*sum;
}
template<std::size_t N> inline BiSourceJet<N> operator/(const BiSourceJet<N>&x,const BiSourceJet<N>&y){return x*bi_inverse(y);}
template<std::size_t N> inline BiSourceJet<N> sqrt(const BiSourceJet<N>&x){
    const Jet x0=x.a[0][0],sx0=sqrt(x0);BiSourceJet<N>h=(x-BiSourceJet<N>{x0})*inverse(x0);
    BiSourceJet<N>sum{Jet(1)},power{Jet(1)};double coefficient=1;
    for(std::size_t k=1;k<=2*N;++k){power=power*h;coefficient*= (1.5-static_cast<double>(k))/static_cast<double>(k);sum=sum+Jet(coefficient)*power;}return sx0*sum;
}
template<std::size_t N> inline BiSourceJet<N> log(const BiSourceJet<N>&x){
    const Jet x0=x.a[0][0];BiSourceJet<N>h=(x-BiSourceJet<N>{x0})*inverse(x0),sum{log(x0)},power{Jet(1)};
    for(std::size_t k=1;k<=2*N;++k){power=power*h;const double coefficient=(k%2?1.0:-1.0)/static_cast<double>(k);sum=sum+Jet(coefficient)*power;}return sum;
}

struct Vec3 { double x, y, z; };
struct JVec3 { Jet x, y, z; };
struct DirectionCosine {
    Jet c{};
    Jet one_minus_c{};
};

inline Vec3 add(Vec3 a, Vec3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 mul(double s, Vec3 a) { return {s*a.x, s*a.y, s*a.z}; }
inline double dot(Vec3 a, Vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Jet dot(const JVec3& a, Vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Jet dot(const JVec3& a, const JVec3& b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline JVec3 div(const JVec3& a, const Jet& s) { return {a.x/s, a.y/s, a.z/s}; }

// Keep the chord-length form of 1-c as a separate Jet.  Constructing only
// c=n.u loses an angular separation below one ulp when the rounded dot
// product is exactly one.  The ordinary dot-product Jet remains preferable
// away from that endpoint; s carries the otherwise unrepresentable separation
// and its observation derivatives for the endpoint branch.
inline DirectionCosine direction_cosine(const JVec3& n, Vec3 u) {
    const Jet dx=n.x-u.x,dy=n.y-u.y,dz=n.z-u.z;
    const Jet s=0.5*(dx*dx+dy*dy+dz*dz);
    return {dot(n,u),s};
}

inline Vec3 radial(double lat, double lon) {
    return {std::cos(lat)*std::cos(lon), std::cos(lat)*std::sin(lon), std::sin(lat)};
}
inline Vec3 north(double lat, double lon) {
    return {-std::sin(lat)*std::cos(lon), -std::sin(lat)*std::sin(lon), std::cos(lat)};
}
inline Vec3 east(double lon) { return {-std::sin(lon), std::cos(lon), 0.0}; }
inline double sinc(double x) {
    const double x2=x*x;
    return std::abs(x)<1.0e-4
        ? 1.0-x2/6.0+x2*x2/120.0-x2*x2*x2/5040.0
        : std::sin(x)/x;
}

inline JVec3 observation_jet(const Point& p) {
    Vec3 n = radial(p.lat, p.lon);
    Vec3 ex = north(p.lat, p.lon);
    Vec3 ey = east(p.lon);
    if(std::abs(std::cos(p.lat))<=64.0*std::numeric_limits<double>::epsilon()){
        const double s=p.lat>=0?1.0:-1.0;
        n={0.0,0.0,s};ex={-s,0.0,0.0};ey={0.0,1.0,0.0};
    }
    const Vec3 ez = n;
    const Jet X = Jet::variable(0.0, 0);
    const Jet Y = Jet::variable(0.0, 1);
    const Jet Z = Jet::variable(0.0, 2);
    return {p.r*n.x + ex.x*X + ey.x*Y + ez.x*Z,
            p.r*n.y + ex.y*X + ey.y*Y + ez.y*Z,
            p.r*n.z + ex.z*X + ey.z*Y + ez.z*Z};
}

inline Jet mean(double a, const Jet& r) {
    if (r.value() < a) return 1.0 - (r*r)/(3.0*a*a);
    return (2.0*a/3.0) / r;
}

inline Jet distance(double a, const Jet& r, const Jet& c) {
    // d_a^2=(a-r)^2+2ar(1-c).  Unlike a^2+r^2-2arc this retains
    // the physical separation when a~=r and c~=1 (Test.tex, Sec. 5.5).
    const Jet ar = Jet(a)-r;
    return sqrt(ar*ar + 2.0*a*r*(1.0-c));
}

inline Jet stable_log_argument(double a,const Jet& r,const Jet& c,
                               const Jet& d) {
    const Jet direct=Jet(a)-r*c+d;
    if(direct.value()>0.25*d.value())return direct;
    // (d+a-rc)(d-a+rc)=r^2(1-c^2).  Select the quotient only when
    // the direct sum loses digits; both forms are algebraically identical.
    return r*r*(1.0-c)*(1.0+c)/(d-Jet(a)+r*c);
}

// Cancellation-free gauge coefficient Q_a in Test.tex, Eq. (Qa_stable_radial).
// It is the vector-potential coefficient that converts the spherical normal
// flux of F0 into its complete-sphere mean.  At r=a the outside formula
// supplies the common one-sided analytic jet; the assembled field has the
// same limit from both sides.
inline Jet q(double a, const Jet& r, const Jet& c) {
    const Jet d = distance(a, r, c);
    if (r.value() < a) {
        const Jet delta=Jet(a)-r;
        const Jet dplus=d+delta;
        return (2.0*a/3.0)*r*r*
               (2.0*a*a*(1.0-c)/dplus-delta)/
               (dplus*(Jet(a)+d+r));
    }
    const Jet delta=r-Jet(a);
    return (2.0*a*a/3.0)*
           (r*d+2.0*delta*(r+a))/
           ((d+delta)*(r+d+a));
}

// Stable logarithm of z2/z1, z_i=r_i-r*c+d_i.  If u_i=r_i-r*c<0,
// z_i is evaluated through z_i=r^2(1-c^2)/(d_i-u_i).  Treating the two
// endpoints separately also covers the mixed-sign interval.
inline Jet stable_log_z_ratio(double r1,double r2,const Jet&r,const Jet&c,
                              const Jet&d1,const Jet&d2) {
    const Jet u1=Jet(r1)-r*c,u2=Jet(r2)-r*c;
    if(u1.value()<0.0) {
        if(u2.value()<0.0)
            return log((d1-u1)/(d2-u2));
        const Jet log_z1=2.0*log(r)+log(1.0-c)+log(1.0+c)-log(d1-u1);
        return log(d2+u2)-log_z1;
    }
    return log((d2+u2)/(d1+u1));
}

inline Jet h(double r1, double r2, const Jet& r, const Jet& c) {
    const Jet d1 = distance(r1, r, c);
    const Jet d2 = distance(r2, r, c);
    return r*(d2-d1)+r*r*c*stable_log_z_ratio(r1,r2,r,c,d1,d2);
}

// Chord-based forms for the sub-ulp neighbourhood of c=1.  In particular,
// these remain finite when r1<r<r2 and a positive angular separation is too
// small for the rounded dot product c=n.u to differ from one.
inline Jet distance_from_s(double a,const Jet&r,const Jet&s) {
    const Jet ar=Jet(a)-r;
    return sqrt(ar*ar+2.0*a*r*s);
}

inline Jet q_from_s(double a,const Jet&r,const Jet&s) {
    const Jet d=distance_from_s(a,r,s);
    if(r.value()<a) {
        const Jet delta=Jet(a)-r,dplus=d+delta;
        return (2.0*a/3.0)*r*r*
               (2.0*a*a*s/dplus-delta)/
               (dplus*(Jet(a)+d+r));
    }
    const Jet delta=r-Jet(a);
    return (2.0*a*a/3.0)*
           (r*d+2.0*delta*(r+a))/
           ((d+delta)*(r+d+a));
}

inline Jet h_from_s(double r1,double r2,const Jet&r,const Jet&s) {
    const Jet d1=distance_from_s(r1,r,s);
    const Jet d2=distance_from_s(r2,r,s);
    const Jet u1=Jet(r1)-r+r*s,u2=Jet(r2)-r+r*s;
    Jet log_ratio;
    if(u1.value()<0.0) {
        if(u2.value()<0.0)
            log_ratio=log((d1-u1)/(d2-u2));
        else {
            const Jet log_z1=2.0*log(r)+log(s)+log(2.0-s)-log(d1-u1);
            log_ratio=log(d2+u2)-log_z1;
        }
    } else {
        log_ratio=log((d2+u2)/(d1+u1));
    }
    return r*(d2-d1)+r*r*(1.0-s)*log_ratio;
}

inline Jet w_hat_from_s(const Prism&prism,const Jet&r,const Jet&s) {
    return (q_from_s(prism.r2,r,s)-q_from_s(prism.r1,r,s)-
            h_from_s(prism.r1,prism.r2,r,s))/(prism.r2-prism.r1);
}

// Tangential side-face coefficient produced directly by the gauge-optimized
// field F_star: W = Q_r2 - Q_r1 - H.  This is algebraically identical to the
// former facewise surface-divergence construction.
inline Jet gauge_boundary_kernel(const Prism& prism, const Jet& r, const Jet& c) {
    return q(prism.r2, r, c) - q(prism.r1, r, c) - h(prism.r1, prism.r2, r, c);
}

// Endpoint primitive used by the cancellation-resistant representation in
// Test.tex, Eq. (cf_U).  The arbitrary logarithmic reference length cancels
// from every source-radius divided difference; one metre is used here.
inline Jet source_primitive(double a, const Jet& r, const Jet& c) {
    const Jet d = distance(a, r, c);
    return r*d + r*r*c*log(stable_log_argument(a,r,c,d));
}

inline Jet u(double a, const Jet& r, const Jet& c) {
    return q(a,r,c)-source_primitive(a,r,c);
}


template<std::size_t N> inline SourceJet<N> source_distance(
    double a,const Jet&r,const SourceJet<N>&c) {
    const SourceJet<N>rs{r},as{Jet(a)};
    return sqrt((as-rs)*(as-rs)+Jet(2.0*a)*rs*(SourceJet<N>{Jet(1.0)}-c));
}

template<std::size_t N> inline SourceJet<N> source_q(
    double a,const Jet&r,const SourceJet<N>&c) {
    const SourceJet<N>rs{r},as{Jet(a)},d=source_distance(a,r,c);
    if(r.value()<a) {
        const SourceJet<N>delta=as-rs,dplus=d+delta;
        return Jet(2.0*a/3.0)*rs*rs*
               (Jet(2.0*a*a)*(SourceJet<N>{Jet(1.0)}-c)/dplus-delta)/
               (dplus*(as+d+rs));
    }
    const SourceJet<N>delta=rs-as;
    return Jet(2.0*a*a/3.0)*
           (rs*d+Jet(2.0)*delta*(rs+as))/
           ((d+delta)*(rs+d+as));
}

template<std::size_t N> inline SourceJet<N> source_primitive(
    double a,const Jet&r,const SourceJet<N>&c) {
    const SourceJet<N>rs{r},as{Jet(a)},d=source_distance(a,r,c);
    SourceJet<N>arg=as-rs*c+d;
    if(arg.a[0].value()<=0.25*d.a[0].value())
        arg=rs*rs*(SourceJet<N>{Jet(1.0)}-c)*
            (SourceJet<N>{Jet(1.0)}+c)/(d-as+rs*c);
    return rs*d+rs*rs*c*log(arg);
}

template<std::size_t N> inline SourceJet<N> source_u(
    double a,const Jet&r,const SourceJet<N>&c) {
    return source_q(a,r,c)-source_primitive(a,r,c);
}

template<std::size_t N> inline SourceJet<N> source_distance_from_s(
    double a,const Jet&r,const SourceJet<N>&s) {
    const SourceJet<N>rs{r},as{Jet(a)};
    return sqrt((as-rs)*(as-rs)+Jet(2.0*a)*rs*s);
}

template<std::size_t N> inline SourceJet<N> source_q_from_s(
    double a,const Jet&r,const SourceJet<N>&s) {
    const SourceJet<N>rs{r},as{Jet(a)},d=source_distance_from_s(a,r,s);
    if(r.value()<a) {
        const SourceJet<N>delta=as-rs,dplus=d+delta;
        return Jet(2.0*a/3.0)*rs*rs*
               (Jet(2.0*a*a)*s/dplus-delta)/
               (dplus*(as+d+rs));
    }
    const SourceJet<N>delta=rs-as;
    return Jet(2.0*a*a/3.0)*
           (rs*d+Jet(2.0)*delta*(rs+as))/
           ((d+delta)*(rs+d+as));
}

template<std::size_t N> inline SourceJet<N> source_primitive_from_s(
    double a,const Jet&r,const SourceJet<N>&s) {
    const SourceJet<N>rs{r},as{Jet(a)},d=source_distance_from_s(a,r,s);
    SourceJet<N>arg=as-rs+rs*s+d;
    if(arg.a[0].value()<=0.25*d.a[0].value())
        arg=rs*rs*s*(SourceJet<N>{Jet(2.0)}-s)/
            (d-as+rs-rs*s);
    return rs*d+rs*rs*(SourceJet<N>{Jet(1.0)}-s)*log(arg);
}

template<std::size_t N> inline SourceJet<N> source_u_from_s(
    double a,const Jet&r,const SourceJet<N>&s) {
    return source_q_from_s(a,r,s)-source_primitive_from_s(a,r,s);
}

inline bool jet_term_small(const Jet&term,const Jet&sum,double factor=128.0){
    const double eps=factor*std::numeric_limits<double>::epsilon();
    for(int i=0;i<=3;++i)for(int j=0;j<=3-i;++j)for(int k=0;k<=3-i-j;++k)
        if(std::abs(term.at(i,j,k))>eps*std::max(1.0,std::abs(sum.at(i,j,k))))return false;
    return true;
}

template<std::size_t N> inline SourceJet<N> radial_mean(
    const SourceJet<N>&a,const Jet&r) {
    const SourceJet<N>rs{r};
    if(r.value()<a.a[0].value())
        return SourceJet<N>{Jet(1.0)}-rs*rs/(Jet(3.0)*a*a);
    return Jet(2.0/3.0)*a/rs;
}

template<std::size_t N> inline SourceJet<N> radial_distance(
    const SourceJet<N>&a,const Jet&r,const Jet&c) {
    const SourceJet<N>rs{r};
    return sqrt((a-rs)*(a-rs)+Jet(2.0)*a*rs*Jet(1.0-c));
}

template<std::size_t N> inline SourceJet<N> radial_q(
    const SourceJet<N>&a,const Jet&r,const Jet&c) {
    const SourceJet<N>rs{r},d=radial_distance(a,r,c);
    if(r.value()<a.a[0].value()) {
        const SourceJet<N>delta=a-rs,dplus=d+delta;
        return Jet(2.0/3.0)*a*rs*rs*
               (Jet(2.0)*a*a*Jet(1.0-c)/dplus-delta)/
               (dplus*(a+d+rs));
    }
    const SourceJet<N>delta=rs-a;
    return Jet(2.0/3.0)*a*a*
           (rs*d+Jet(2.0)*delta*(rs+a))/
           ((d+delta)*(rs+d+a));
}

template<std::size_t N> inline SourceJet<N> radial_primitive(
    const SourceJet<N>&a,const Jet&r,const Jet&c) {
    const SourceJet<N>rs{r},d=radial_distance(a,r,c);
    SourceJet<N>arg=a-rs*Jet(c)+d;
    if(arg.a[0].value()<=0.25*d.a[0].value())
        arg=rs*rs*Jet(1.0-c)*Jet(1.0+c)/(d-a+rs*Jet(c));
    return rs*d+rs*rs*Jet(c)*log(arg);
}

template<std::size_t N> inline SourceJet<N> radial_u(
    const SourceJet<N>&a,const Jet&r,const Jet&c) {
    return radial_q(a,r,c)-radial_primitive(a,r,c);
}

template<std::size_t N> inline BiSourceJet<N> bi_distance(
    const BiSourceJet<N>&a,const Jet&r,const BiSourceJet<N>&c){
    const BiSourceJet<N>rs{r};return sqrt((a-rs)*(a-rs)+Jet(2)*a*rs*(BiSourceJet<N>{Jet(1)}-c));
}
template<std::size_t N> inline BiSourceJet<N> bi_q(
    const BiSourceJet<N>&a,const Jet&r,const BiSourceJet<N>&c){
    const BiSourceJet<N>rs{r},d=bi_distance(a,r,c);
    if(r.value()<a.a[0][0].value()){
        const BiSourceJet<N>delta=a-rs,dplus=d+delta;
        return Jet(2.0/3.0)*a*rs*rs*
               (Jet(2.0)*a*a*(BiSourceJet<N>{Jet(1.0)}-c)/dplus-delta)/
               (dplus*(a+d+rs));
    }
    const BiSourceJet<N>delta=rs-a;
    return Jet(2.0/3.0)*a*a*
           (rs*d+Jet(2.0)*delta*(rs+a))/
           ((d+delta)*(rs+d+a));
}
template<std::size_t N> inline BiSourceJet<N> bi_primitive(
    const BiSourceJet<N>&a,const Jet&r,const BiSourceJet<N>&c){
    const BiSourceJet<N>rs{r},d=bi_distance(a,r,c);BiSourceJet<N>arg=a-rs*c+d;
    if(arg.a[0][0].value()<=0.25*d.a[0][0].value())
        arg=rs*rs*(BiSourceJet<N>{Jet(1)}-c)*
            (BiSourceJet<N>{Jet(1)}+c)/(d-a+rs*c);
    return rs*d+rs*rs*c*log(arg);
}
template<std::size_t N> inline BiSourceJet<N> bi_u(
    const BiSourceJet<N>&a,const Jet&r,const BiSourceJet<N>&c){return bi_q(a,r,c)-bi_primitive(a,r,c);}

inline Jet double_divided_difference_ac(const Prism&prism,const Jet&r,
                                        const Jet&cm,const Jet&cp){
    const Jet cmid=0.5*(cm+cp),hc=0.5*(cp-cm);
    auto local=[&](double lo,double hi){
        const double amid=0.5*(lo+hi),ha=0.5*(hi-lo);
        const auto aa=BiSourceJet<3>::variable_a(Jet(amid));
        const auto cc=BiSourceJet<3>::variable_c(cmid);const auto series=bi_u(aa,r,cc);
        Jet sum;double hap=1.0;
        for(std::size_t i=1;i<=3;i+=2){Jet hcp=Jet(1);for(std::size_t j=1;j<=3;j+=2){sum=sum+hap*hcp*series.a[i][j];hcp=hcp*hc*hc;}hap*=ha*ha;}
        return sum;
    };
    const double rv=r.value(),width=prism.r2-prism.r1;
    if(prism.r1<rv&&rv<prism.r2)
        return ((rv-prism.r1)*local(prism.r1,rv)+(prism.r2-rv)*local(rv,prism.r2))/width;
    return local(prism.r1,prism.r2);
}

inline Jet radial_u_divided_difference(double a1,double a2,
                                       const Jet&r,const Jet&c) {
    const double width=a2-a1,scale=std::max({1.0,std::abs(a1),std::abs(a2)});
    if(std::abs(width)>1.0e-4*scale)
        return (u(a2,r,c)-u(a1,r,c))/width;
    auto local=[&](double lo,double hi){
        const double midpoint=0.5*(lo+hi),half=0.5*(hi-lo),h2=half*half;
        const SourceJet<9>as=SourceJet<9>::variable(Jet(midpoint));
        const SourceJet<9>series=radial_u(as,r,c);
        Jet sum=series.a[1];double hp=h2;
        for(std::size_t order=3;order<=9;order+=2){const Jet term=hp*series.a[order];sum=sum+term;if(jet_term_small(term,sum)&&order>=5)break;hp*=h2;}
        return sum;
    };
    const double rv=r.value();
    if(a1<rv&&rv<a2)
        return ((rv-a1)*local(a1,rv)+(a2-rv)*local(rv,a2))/width;
    return local(a1,a2);
}

inline Jet radial_f_divided_difference(double a1,double a2,
                                       const Jet&r) {
    const double width=a2-a1,scale=std::max({1.0,std::abs(a1),std::abs(a2)});
    if(std::abs(width)>1.0e-4*scale)
        return (a2*a2*mean(a2,r)-a1*a1*mean(a1,r))/width;
    auto local=[&](double lo,double hi){
        const double midpoint=0.5*(lo+hi),half=0.5*(hi-lo),h2=half*half;
        const SourceJet<9>as=SourceJet<9>::variable(Jet(midpoint));
        const SourceJet<9>series=as*as*radial_mean(as,r);
        Jet sum=series.a[1];double hp=h2;
        for(std::size_t order=3;order<=9;order+=2){const Jet term=hp*series.a[order];sum=sum+term;if(jet_term_small(term,sum)&&order>=5)break;hp*=h2;}
        return sum;
    };
    const double rv=r.value();
    if(a1<rv&&rv<a2)
        return ((rv-a1)*local(a1,rv)+(a2-rv)*local(rv,a2))/width;
    return local(a1,a2);
}

inline Jet w_hat(const Prism& prism,const Jet&r,const Jet&c);

inline double fact(int n){double v=1;for(int i=2;i<=n;++i)v*=i;return v;}
inline double power_log(double v,int k){
    return k==-1?std::log(v):std::pow(v,static_cast<double>(k+1))/static_cast<double>(k+1);
}

inline double physical_moment_collinear(const Prism& prism,double r,
                                              int endpoint,int n,int p){
    double sum1=0,sum2=0;
    auto choose=[](int nn,int kk){double v=1;for(int i=1;i<=kk;++i)v*=static_cast<double>(nn+1-i)/i;return v;};
    if(endpoint>0){
        if(prism.r1<=r&&r<=prism.r2)return std::numeric_limits<double>::infinity();
        const int sigma=prism.r1>r?1:-1;
        const double v1=std::abs(prism.r1-r),v2=std::abs(prism.r2-r);
        for(int j=0;j<=p;++j){
            const double coefficient=sigma*choose(p,j)*std::pow(r,p-j)*
                (j%2&&sigma<0?-1.0:1.0);
            sum1+=coefficient*power_log(v1,j-n);sum2+=coefficient*power_log(v2,j-n);
        }
    }else{
        const double v1=prism.r1+r,v2=prism.r2+r;
        for(int j=0;j<=p;++j){
            const double coefficient=choose(p,j)*std::pow(-r,p-j);
            sum1+=coefficient*power_log(v1,j-n);sum2+=coefficient*power_log(v2,j-n);
        }
    }
    return sum2-sum1;
}

inline double h_collinear(const Prism& prism,double r,int endpoint){
    const double dr=prism.r2-prism.r1;
    if(endpoint<0)return r*(dr-r*std::log1p(dr/(r+prism.r1)));
    if(r<prism.r1)return r*(dr+r*std::log1p(dr/(prism.r1-r)));
    if(r>prism.r2)return r*(r*std::log1p(dr/(r-prism.r2))-dr);
    return std::numeric_limits<double>::infinity();
}

inline std::array<std::array<double,4>,4> h_endpoint_partials(
    const Prism& prism,double r,int endpoint){
    std::array<std::array<double,4>,4> h{};const double c=endpoint;
    auto m=[&](int n,int p){return physical_moment_collinear(prism,r,endpoint,n,p);};
    h[0][0]=h_collinear(prism,r,endpoint);
    h[1][0]=m(1,1)-r*r*m(3,1)+r*c*m(3,2);
    h[0][1]=r*r*m(3,2);
    h[2][0]=-3*r*m(3,1)+2*c*m(3,2)+3*r*r*r*m(5,1)-6*c*r*r*m(5,2)+3*c*c*r*m(5,3);
    h[1][1]=2*r*m(3,2)-3*r*r*r*m(5,2)+3*c*r*r*m(5,3);
    h[0][2]=3*r*r*r*m(5,3);
    h[0][3]=15*r*r*r*r*m(7,4);
    h[1][2]=9*r*r*m(5,3)-15*r*r*r*r*m(7,3)+15*c*r*r*r*m(7,4);
    h[2][1]=2*m(3,2)-15*r*r*m(5,2)+12*c*r*m(5,3)+15*r*r*r*r*m(7,2)
             -30*c*r*r*r*m(7,3)+15*c*c*r*r*m(7,4);
    h[3][0]=-3*m(3,1)+18*r*r*m(5,1)-27*c*r*m(5,2)+9*c*c*m(5,3)
             -15*r*r*r*r*m(7,1)+45*c*r*r*r*m(7,2)-45*c*c*r*r*m(7,3)+15*c*c*c*r*m(7,4);
    return h;
}

inline std::array<std::array<double,4>,4> q_endpoint_partials(
    double a,double r0,int endpoint){
    // The endpoint Taylor table is an auxiliary expansion in (r,c), not in
    // the requested observation derivative alone.  It needs every Jet
    // coefficient through total degree three even when the caller is using a
    // component-selective sparse context.  Restore the caller's sparse mask
    // automatically when this local table has been formed.
    const ScopedJetContext endpoint_context(20);
    std::array<std::array<double,4>,4> out{};
    const Jet rr=Jet::variable(r0,0),cc=Jet::variable(static_cast<double>(endpoint),1);
    const Jet d=distance(a,rr,cc);
    const Jet f=(Jet(a)-rr*cc)/d;
    const Jet g=f-mean(a,rr);
    for(int m=0;m<=3;++m)for(int n=0;n<=3-m;++n){
        double sum=0;
        for(int k=0;k<=n;++k){
            const double gmk=g.at(m,k,0)*fact(m)*fact(k);
            const double sign=(endpoint>0&&k%2)?-1.0:1.0;
            sum+=sign*gmk/(fact(k+1)*std::pow(2.0,n-k+1));
        }
        double pref=a*a*fact(n);
        if(endpoint>0&&n%2==0)pref=-pref;
        out[m][n]=pref*sum;
    }
    return out;
}

inline Jet w_hat_endpoint_delta(const Prism& prism,const Jet&r,
                                const Jet&dc,int endpoint){
    const double r0=r.value();
    const auto q1=q_endpoint_partials(prism.r1,r0,endpoint);
    const auto q2=q_endpoint_partials(prism.r2,r0,endpoint);
    const auto hh=h_endpoint_partials(prism,r0,endpoint);
    const Jet dr=r-r0;
    std::array<Jet,4> rp{Jet(1),dr,dr*dr,dr*dr*dr};
    std::array<Jet,4> cp{Jet(1),dc,dc*dc,dc*dc*dc};
    Jet result;
    for(int m=0;m<=3;++m)for(int n=0;n<=3-m;++n)
        result=result+(q2[m][n]-q1[m][n]-hh[m][n])/(fact(m)*fact(n))*rp[m]*cp[n];
    return result/(prism.r2-prism.r1);
}

inline Jet w_hat_endpoint(const Prism& prism,const Jet&r,const Jet&c,int endpoint){
    return w_hat_endpoint_delta(
        prism,r,c-static_cast<double>(endpoint),endpoint);
}

// Divided difference of the endpoint Taylor polynomial.  The polynomial
// identity (x^n-y^n)/(x-y)=sum x^(n-1-k)y^k avoids subtracting two almost
// identical endpoint kernels.
inline Jet w_hat_endpoint_divided_difference_c(
    const Prism& prism,const Jet&r,const Jet&dcm,const Jet&dcp,int endpoint){
    const double r0=r.value();
    const auto q1=q_endpoint_partials(prism.r1,r0,endpoint);
    const auto q2=q_endpoint_partials(prism.r2,r0,endpoint);
    const auto hh=h_endpoint_partials(prism,r0,endpoint);
    const Jet dr=r-r0;
    std::array<Jet,4> rp{Jet(1),dr,dr*dr,dr*dr*dr};
    std::array<Jet,4> pm{Jet(1),dcm,dcm*dcm,dcm*dcm*dcm};
    std::array<Jet,4> pp{Jet(1),dcp,dcp*dcp,dcp*dcp*dcp};
    Jet result;
    for(int m=0;m<=3;++m)for(int n=1;n<=3-m;++n){
        Jet quotient;
        for(int k=0;k<n;++k) quotient=quotient+pp[n-1-k]*pm[k];
        result=result+(q2[m][n]-q1[m][n]-hh[m][n])/
                      (fact(m)*fact(n))*rp[m]*quotient;
    }
    return result/(prism.r2-prism.r1);
}

// Local c-divided difference of W_hat.  Coefficients 1,3,5 are the
// normalized odd derivatives in Eq. (cf_divided_difference_finite).
inline Jet w_hat_divided_difference_c(const Prism& prism,const Jet&r,
                                      const Jet&cm,const Jet&cp) {
    const Jet delta=cp-cm;
    const double radial_scale=std::max({1.0,std::abs(prism.r1),std::abs(prism.r2)});
    if(std::abs(prism.r2-prism.r1)<=1.0e-8*radial_scale&&
       std::abs(delta.value())<=1.0e-6*std::max(1.0,std::abs(0.5*(cp.value()+cm.value()))))
        return double_divided_difference_ac(prism,r,cm,cp);
    if(std::abs(delta.value())>1.0e-4*std::max(1.0,std::abs(0.5*(cp.value()+cm.value()))))
        return (w_hat(prism,r,cp)-w_hat(prism,r,cm))/delta;
    if(std::abs(delta.value())<=32.0*std::numeric_limits<double>::epsilon()){
        constexpr std::size_t Z=3;const Jet midpoint=0.5*(cp+cm),half=0.5*delta;
        const SourceJet<Z>cs=SourceJet<Z>::variable(midpoint);
        const SourceJet<Z>series=(source_u(prism.r2,r,cs)-source_u(prism.r1,r,cs))/
                                  SourceJet<Z>{Jet(prism.r2-prism.r1)};
        return series.a[1]+half*half*series.a[3];
    }
    constexpr std::size_t N=9;
    const Jet midpoint=0.5*(cp+cm),half=0.5*delta;
    const SourceJet<N>cs=SourceJet<N>::variable(midpoint);
    const SourceJet<N>series=(source_u(prism.r2,r,cs)-source_u(prism.r1,r,cs)) /
                              SourceJet<N>{Jet(prism.r2-prism.r1)};
    const Jet h2=half*half;Jet sum=series.a[1],hp=h2;
    for(std::size_t order=3;order<=N;order+=2){const Jet term=hp*series.a[order];sum=sum+term;if(jet_term_small(term,sum)&&order>=5)break;hp=hp*h2;}
    return sum;
}

// The gauge boundary coefficient satisfies W = dr * W_hat.  Keeping dr
// outside all angular differences is the first of the three exact
// cancellation removals in Sec. 5.5 of the manuscript.
inline Jet w_hat(const Prism& prism, const Jet& r, const Jet& c) {
    const double endpoint_distance=1.0-std::abs(c.value());
    if(std::abs(endpoint_distance)<=32.0*std::numeric_limits<double>::epsilon()&&
       r.value()!=prism.r1&&r.value()!=prism.r2&&
       !(c.value()>0&&prism.r1<=r.value()&&r.value()<=prism.r2))
        return w_hat_endpoint(prism,r,c,c.value()>=0?1:-1);
    return radial_u_divided_difference(prism.r1,prism.r2,r,c);
}

inline bool positive_endpoint_regular(const Prism& prism,const Jet&r) {
    return r.value()!=prism.r1&&r.value()!=prism.r2&&
           !(prism.r1<=r.value()&&r.value()<=prism.r2);
}

// Chord-aware W_hat.  In the normal range the established radial divided
// difference is used unchanged.  Only the sub-ulp neighbourhood of c=1
// switches to the endpoint Taylor table and uses dc=-s directly.
inline Jet w_hat(const Prism& prism,const Jet&r,
                 const DirectionCosine& angle) {
    if(angle.one_minus_c.value()>=0.0&&
       angle.one_minus_c.value()<=32.0*std::numeric_limits<double>::epsilon()) {
        if(positive_endpoint_regular(prism,r))
            return w_hat_endpoint_delta(prism,r,-angle.one_minus_c,1);
        return w_hat_from_s(prism,r,angle.one_minus_c);
    }
    return w_hat(prism,r,angle.c);
}

inline Jet w_hat_divided_difference_c_from_s(
    const Prism&prism,const Jet&r,const Jet&sm,const Jet&sp) {
    const Jet delta_c=sm-sp;
    const double scale=std::max(1.0,std::abs(0.5*(sm.value()+sp.value())));
    if(std::abs(delta_c.value())>1.0e-4*scale)
        return (w_hat_from_s(prism,r,sp)-w_hat_from_s(prism,r,sm))/delta_c;

    constexpr std::size_t N=9;
    const Jet midpoint=0.5*(sp+sm),half=0.5*(sp-sm);
    const SourceJet<N>ss=SourceJet<N>::variable(midpoint);
    const SourceJet<N>series=
        (source_u_from_s(prism.r2,r,ss)-source_u_from_s(prism.r1,r,ss))/
        SourceJet<N>{Jet(prism.r2-prism.r1)};
    const Jet h2=half*half;
    Jet sum=-series.a[1],hp=h2;
    for(std::size_t order=3;order<=N;order+=2) {
        const Jet term=-hp*series.a[order];
        sum=sum+term;
        if(jet_term_small(term,sum)&&order>=5)break;
        hp=hp*h2;
    }
    return sum;
}

inline Jet w_hat_divided_difference_c(
    const Prism& prism,const Jet&r,const DirectionCosine& am,
    const DirectionCosine& ap) {
    const Jet delta=ap.c-am.c;
    const double scale=std::max(
        1.0,std::abs(0.5*(ap.c.value()+am.c.value())));
    if(std::abs(delta.value())>1.0e-4*scale)
        return (w_hat(prism,r,ap)-w_hat(prism,r,am))/delta;
    if(am.one_minus_c.value()>=0.0&&ap.one_minus_c.value()>=0.0&&
       std::max(am.one_minus_c.value(),ap.one_minus_c.value())<=
           32.0*std::numeric_limits<double>::epsilon()) {
        if(positive_endpoint_regular(prism,r))
            return w_hat_endpoint_divided_difference_c(
                prism,r,-am.one_minus_c,-ap.one_minus_c,1);
        return w_hat_divided_difference_c_from_s(
            prism,r,am.one_minus_c,ap.one_minus_c);
    }
    return w_hat_divided_difference_c(prism,r,am.c,ap.c);
}

inline std::array<int, 3> multi_index(Component c) {
    switch (c) {
    case Component::V: return {0,0,0};
    case Component::Vx: return {1,0,0}; case Component::Vy: return {0,1,0}; case Component::Vz: return {0,0,1};
    case Component::Vxx: return {2,0,0}; case Component::Vxy: return {1,1,0}; case Component::Vxz: return {1,0,1};
    case Component::Vyy: return {0,2,0}; case Component::Vyz: return {0,1,1}; case Component::Vzz: return {0,0,2};
    case Component::Vxxx: return {3,0,0}; case Component::Vxxy: return {2,1,0}; case Component::Vxxz: return {2,0,1};
    case Component::Vxyy: return {1,2,0}; case Component::Vxyz: return {1,1,1}; case Component::Vxzz: return {1,0,2};
    case Component::Vyyy: return {0,3,0}; case Component::Vyyz: return {0,2,1}; case Component::Vyzz: return {0,1,2};
    case Component::Vzzz: return {0,0,3};
    }
    return {0,0,0};
}

inline double factorial(int n) { return n < 2 ? 1.0 : (n == 2 ? 2.0 : 6.0); }
inline double derivative(const Jet& x, Component c) {
    const auto m = multi_index(c);
    return x.at(m[0],m[1],m[2])*
           factorial(m[0])*factorial(m[1])*factorial(m[2]);
}


using Components = std::array<double,20>;

inline Components derivatives(const Jet& x) {
    Components out{};
    for (int i=0;i<20;++i)
        out[static_cast<std::size_t>(i)]=derivative(x,static_cast<Component>(i));
    return out;
}

inline void source_geometry(const Edge& edge, double parameter,
                            Vec3& u, Vec3& tangent, double& alpha) {
    if (edge.kind == EdgeKind::Latitude) {
        u = radial(edge.fixed, parameter);
        tangent = north(edge.fixed, parameter);
        alpha = edge.sign * std::cos(edge.fixed);
    } else {
        u = radial(parameter, edge.fixed);
        tangent = east(edge.fixed);
        alpha = static_cast<double>(edge.sign);
    }
}

} // namespace detail

inline std::array<Edge,4> edges(const Prism& p) {
    return {Edge{EdgeKind::Latitude,p.lat1,p.lon1,p.lon2,-1},
            Edge{EdgeKind::Latitude,p.lat2,p.lon1,p.lon2,+1},
            Edge{EdgeKind::Longitude,p.lon1,p.lat1,p.lat2,-1},
            Edge{EdgeKind::Longitude,p.lon2,p.lat1,p.lat2,+1}};
}

// One edge integrand for exactly one requested physical component.
inline double kernel(Component component, const Prism& prism,
                          const Point& point, const Edge& edge,
                          double parameter) {
    using namespace detail;
    const ScopedJetContext context(static_cast<std::size_t>(component));
    const JVec3 x = observation_jet(point);
    const Jet r = sqrt(dot(x,x));
    const JVec3 n = div(x,r);
    Vec3 u{}, tangent{};
    double alpha = 0.0;
    source_geometry(edge, parameter, u, tangent, alpha);
    const DirectionCosine c = direction_cosine(n,u);
    const Jet b = alpha * dot(n,tangent);
    return derivative(b*(prism.r2-prism.r1)*w_hat(prism,r,c), component);
}

inline detail::Components kernel_all(const Prism& prism,const Point& point,
                                     const Edge& edge,double parameter) {
    using namespace detail;
    const ScopedJetContext context(20);
    const JVec3 x=observation_jet(point);const Jet r=sqrt(dot(x,x));const JVec3 n=div(x,r);
    Vec3 u0{},t{};double alpha=0;source_geometry(edge,parameter,u0,t,alpha);
    const DirectionCosine c=direction_cosine(n,u0);
    const Jet b=alpha*dot(n,t);
    return derivatives(b*(prism.r2-prism.r1)*w_hat(prism,r,c));
}

inline double explicit_term(Component component, const Prism& prism,
                                 const Point& point) {
    using namespace detail;
    const ScopedJetContext context(static_cast<std::size_t>(component));
    const JVec3 x = observation_jet(point);
    const Jet r = sqrt(dot(x,x));
    const double omega = (prism.lon2-prism.lon1) *
                              (std::sin(prism.lat2)-std::sin(prism.lat1));
    // Constant radial flux of the gauge-optimized field: C0.
    const Jet c0 = omega * (prism.r2*prism.r2*mean(prism.r2,r)
                           -prism.r1*prism.r1*mean(prism.r1,r));
    return derivative(c0, component);
}

// The two opposite edges are evaluated as an average times a divided
// difference (Test.tex, Eqs. cf_product_difference, cf_Kphi and cf_Klambda).
// This routine returns the *scaled* paired kernel, i.e. the angular cell width
// and radial thickness have already been divided out.
namespace detail {
struct PreparedPairedKernel {
    Jet r{};
    JVec3 n{};
    EdgeKind kind=EdgeKind::Latitude;
    double width=0.0;
    double lower=0.0,upper=0.0,mid=0.0,half_width=0.0;
    double sin_lower=0.0,cos_lower=1.0;
    double sin_upper=0.0,cos_upper=1.0;
    double sin_mid=0.0,cos_mid=1.0;
    double sinc_half=1.0,sinc_full=1.0;
    Vec3 tangent_lower{},tangent_upper{},horizontal_mid{},east_mid{};
};
} // namespace detail

inline detail::PreparedPairedKernel prepare_paired_kernel(
    const Prism& prism,const Point& point,EdgeKind kind){
    using namespace detail;
    PreparedPairedKernel prepared;
    const JVec3 x=observation_jet(point);
    prepared.r=sqrt(dot(x,x));
    prepared.n=div(x,prepared.r);
    prepared.kind=kind;
    prepared.lower=kind==EdgeKind::Latitude?prism.lat1:prism.lon1;
    prepared.upper=kind==EdgeKind::Latitude?prism.lat2:prism.lon2;
    prepared.width=prepared.upper-prepared.lower;
    prepared.mid=0.5*(prepared.lower+prepared.upper);
    prepared.half_width=0.5*prepared.width;
    prepared.sin_lower=std::sin(prepared.lower);
    prepared.cos_lower=std::cos(prepared.lower);
    prepared.sin_upper=std::sin(prepared.upper);
    prepared.cos_upper=std::cos(prepared.upper);
    prepared.sin_mid=std::sin(prepared.mid);
    prepared.cos_mid=std::cos(prepared.mid);
    prepared.sinc_half=sinc(prepared.half_width);
    prepared.sinc_full=sinc(prepared.width);
    if(kind==EdgeKind::Longitude){
        prepared.tangent_lower=east(prepared.lower);
        prepared.tangent_upper=east(prepared.upper);
        prepared.horizontal_mid={-prepared.cos_mid,-prepared.sin_mid,0.0};
        prepared.east_mid=east(prepared.mid);
    }
    return prepared;
}

inline detail::Jet paired_kernel_jet(const Prism& prism,
                                     const detail::PreparedPairedKernel& prepared,
                                     double parameter) {
    using namespace detail;
    const Jet& r=prepared.r;
    const JVec3& n=prepared.n;
    const EdgeKind kind=prepared.kind;

    Vec3 um{},up{},tm{},tp{};
    if(kind==EdgeKind::Latitude) {
        const double sp=std::sin(parameter),cp=std::cos(parameter);
        um={prepared.cos_lower*cp,prepared.cos_lower*sp,prepared.sin_lower};
        up={prepared.cos_upper*cp,prepared.cos_upper*sp,prepared.sin_upper};
        tm={-prepared.sin_lower*cp,-prepared.sin_lower*sp,prepared.cos_lower};
        tp={-prepared.sin_upper*cp,-prepared.sin_upper*sp,prepared.cos_upper};
    } else {
        const double sp=std::sin(parameter),cp=std::cos(parameter);
        um={cp*prepared.cos_lower,cp*prepared.sin_lower,sp};
        up={cp*prepared.cos_upper,cp*prepared.sin_upper,sp};
        tm=prepared.tangent_lower;
        tp=prepared.tangent_upper;
    }
    const DirectionCosine cm=direction_cosine(n,um);
    const DirectionCosine cp=direction_cosine(n,up);
    const Jet wm=w_hat(prism,r,cm), wp=w_hat(prism,r,cp);
    Jet am,ap;
    if(kind==EdgeKind::Latitude) {
        am=std::cos(prism.lat1)*dot(n,tm);
        ap=std::cos(prism.lat2)*dot(n,tp);
    } else {
        am=dot(n,tm);
        ap=dot(n,tp);
    }
    Jet dd_c,dd_a;
    if(kind==EdgeKind::Latitude) {
        const double sp=std::sin(parameter),cos_parameter=std::cos(parameter);
        const Vec3 radial_mid{prepared.cos_mid*cos_parameter,prepared.cos_mid*sp,prepared.sin_mid};
        const Vec3 north_mid{-prepared.sin_mid*cos_parameter,-prepared.sin_mid*sp,prepared.cos_mid};
        const Jet c_mid=dot(n,radial_mid);
        const Jet c_prime=dot(n,north_mid);
        dd_c=prepared.sinc_half*c_prime;
        dd_a=prepared.sinc_full*(-prepared.cos_mid*c_mid-prepared.sin_mid*c_prime);
    } else {
        dd_a=prepared.sinc_half*dot(n,prepared.horizontal_mid);
        dd_c=prepared.sinc_half*std::cos(parameter)*dot(n,prepared.east_mid);
    }
    const Jet average_a=0.5*(ap+am);
    const Jet average_w=0.5*(wp+wm);
    const Jet dd_w=dd_c*w_hat_divided_difference_c(prism,r,cm,cp);
    return average_a*dd_w+average_w*dd_a;
}

inline detail::Jet paired_kernel_jet(const Prism& prism,const Point& point,
                                     EdgeKind kind,double parameter){
    return paired_kernel_jet(prism,prepare_paired_kernel(prism,point,kind),parameter);
}

inline double paired_kernel_prepared(Component component,const Prism& prism,
                                     const detail::PreparedPairedKernel& prepared,
                                     double parameter){
    const detail::ScopedJetContext context(static_cast<std::size_t>(component));
    return detail::derivative(paired_kernel_jet(prism,prepared,parameter),component);
}

inline detail::Components paired_kernel_all(const Prism& prism, const Point& point,
                                            EdgeKind kind, double parameter) {
    const detail::ScopedJetContext context(20);
    return detail::derivatives(paired_kernel_jet(prism,point,kind,parameter));
}

// Component-selective form used by the independent adaptive integrators.  It
// evaluates the same cancellation-resistant paired Jet as paired_kernel_all,
// but extracts only the requested physical derivative.
inline double paired_kernel(Component component,const Prism& prism,const Point& point,
                            EdgeKind kind,double parameter) {
    const detail::ScopedJetContext context(static_cast<std::size_t>(component));
    return detail::derivative(paired_kernel_jet(prism,point,kind,parameter),component);
}

inline detail::Components scaled_explicit_all(const Prism& prism,const Point& point) {
    using namespace detail;
    const ScopedJetContext context(20);
    const JVec3 x=observation_jet(point);
    const Jet r=sqrt(dot(x,x));
    const double hc=0.5*(prism.lat2-prism.lat1);
    const double latc=0.5*(prism.lat1+prism.lat2);
    const double sinc=detail::sinc(hc);
    const Jet fdd=radial_f_divided_difference(prism.r1,prism.r2,r);
    return derivatives(fdd*(std::cos(latc)*sinc));
}

inline double scaled_explicit(Component component,const Prism& prism,const Point& point) {
    using namespace detail;
    const ScopedJetContext context(static_cast<std::size_t>(component));
    const JVec3 x=observation_jet(point);
    const Jet r=sqrt(dot(x,x));
    const double hc=0.5*(prism.lat2-prism.lat1);
    const double latc=0.5*(prism.lat1+prism.lat2);
    const double sinc=detail::sinc(hc);
    const Jet fdd=radial_f_divided_difference(prism.r1,prism.r2,r);
    return derivative(fdd*(std::cos(latc)*sinc),component);
}

} // namespace Tesseroid
