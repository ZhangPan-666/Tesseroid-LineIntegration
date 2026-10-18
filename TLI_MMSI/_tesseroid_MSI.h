#pragma once

#include "_TLI_estimate.h"
#include "_integral2.h"
#include "_tesseroid_surfacekernel.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Tesseroid {

enum TesseroidMSIComponent : std::uint32_t {
    MSI_V = 1u << 0,
    MSI_Vx = 1u << 1,
    MSI_Vy = 1u << 2,
    MSI_Vz = 1u << 3,
    MSI_Vxx = 1u << 4,
    MSI_Vxy = 1u << 5,
    MSI_Vyy = 1u << 6,
    MSI_Vzx = 1u << 7,
    MSI_Vzy = 1u << 8,
    MSI_Vzz = 1u << 9,
    MSI_Vxxx = 1u << 10,
    MSI_Vxxy = 1u << 11,
    MSI_Vxxz = 1u << 12,
    MSI_Vxyz = 1u << 13,
    MSI_Vyyx = 1u << 14,
    MSI_Vyyy = 1u << 15,
    MSI_Vyyz = 1u << 16,
    MSI_Vzzx = 1u << 17,
    MSI_Vzzy = 1u << 18,
    MSI_Vzzz = 1u << 19,
    MSI_All = (1u << 20) - 1u
};

namespace tli_mmsi_detail {

constexpr double degree = PI / 180.0;
using Values = std::array<double, 20>;

template <class T>
inline T clamp_value(T value, T lower, T upper) {
    return value < lower ? lower : (value > upper ? upper : value);
}

inline double smoothstep5(double t) {
    t = clamp_value(t, 0.0, 1.0);
    return t * t * t * (10.0 + t * (-15.0 + 6.0 * t));
}

inline thread_local std::size_t requested_component_context = 20;
inline thread_local std::uint32_t requested_component_mask =
    (1u << 20) - 1u;

class ScopedComponentContext {
  public:
    explicit ScopedComponentContext(
        std::size_t component,
        std::uint32_t component_mask = (1u << 20) - 1u)
        : previous_(requested_component_context),
          previous_mask_(requested_component_mask) {
        requested_component_context = component;
        requested_component_mask = component_mask;
    }
    ~ScopedComponentContext() {
        requested_component_context = previous_;
        requested_component_mask = previous_mask_;
    }

  private:
    std::size_t previous_;
    std::uint32_t previous_mask_;
};

struct Grid {
    double lon_first_deg = 0.0;
    double lon_last_deg = 0.0;
    double lat_first_deg = 0.0;
    double lat_last_deg = 0.0;
    double dlon_deg = 0.0;
    double dlat_deg = 0.0;
    std::size_t nlon = 0;
    std::size_t nlat = 0;
    bool global_lon = false;

    double lon_west_deg() const { return lon_first_deg - 0.5 * dlon_deg; }
    double lon_east_deg() const { return lon_last_deg + 0.5 * dlon_deg; }
    double lat_south_deg() const {
        return std::max(-90.0, lat_first_deg - 0.5 * dlat_deg);
    }
    double lat_north_deg() const {
        return std::min(90.0, lat_last_deg + 0.5 * dlat_deg);
    }
};

inline Grid make_grid(double lons, double lone, double lats, double late,
                      double dlon, double dlat, std::size_t model_size) {
    if (!(dlon > 0.0 && dlat > 0.0 && lone >= lons && late >= lats))
        throw std::invalid_argument("TLI-MMSI requires positive grid spacing and ordered bounds");

    Grid g;
    g.lon_first_deg = lons;
    g.lon_last_deg = lone;
    g.lat_first_deg = lats;
    g.lat_last_deg = late;
    g.dlon_deg = dlon;
    g.dlat_deg = dlat;
    g.nlon = static_cast<std::size_t>(std::llround((lone - lons) / dlon)) + 1;
    g.nlat = static_cast<std::size_t>(std::llround((late - lats) / dlat)) + 1;
    if (g.nlon == 0 || g.nlat == 0 || g.nlon * g.nlat != model_size)
        throw std::invalid_argument("TLI-MMSI grid dimensions do not match source array size");

    const double span = g.lon_east_deg() - g.lon_west_deg();
    g.global_lon = std::abs(span - 360.0) <=
                   128.0 * std::numeric_limits<double>::epsilon() *
                       std::max(360.0, std::abs(span));
    return g;
}

inline std::size_t index(const Grid& g, std::size_t ilon, std::size_t ilat) {
    return ilon * g.nlat + ilat;
}

inline double interpolate_row(const Grid& g, const std::vector<double>& f,
                              std::size_t row, double lon_deg) {
    double u = (lon_deg - g.lon_first_deg) / g.dlon_deg;
    if (g.global_lon) {
        u = std::fmod(u, static_cast<double>(g.nlon));
        if (u < 0.0)
            u += static_cast<double>(g.nlon);
        const std::size_t i0 = static_cast<std::size_t>(std::floor(u)) % g.nlon;
        const std::size_t i1 = (i0 + 1) % g.nlon;
        const double t = u - std::floor(u);
        return std::fma(t, f[index(g, i1, row)] - f[index(g, i0, row)],
                        f[index(g, i0, row)]);
    }

    if (u <= 0.0)
        return f[index(g, 0, row)];
    if (u >= static_cast<double>(g.nlon - 1))
        return f[index(g, g.nlon - 1, row)];
    const std::size_t i0 = static_cast<std::size_t>(std::floor(u));
    const double t = u - static_cast<double>(i0);
    return std::fma(t, f[index(g, i0 + 1, row)] - f[index(g, i0, row)],
                    f[index(g, i0, row)]);
}

inline double interpolate_column(const Grid& g, const std::vector<double>& f,
                                 std::size_t column, double lat_deg) {
    double u = (lat_deg - g.lat_first_deg) / g.dlat_deg;
    if (u <= 0.0)
        return f[index(g, column, 0)];
    if (u >= static_cast<double>(g.nlat - 1))
        return f[index(g, column, g.nlat - 1)];

    const std::size_t j0 = static_cast<std::size_t>(std::floor(u));
    const double lat0 = (g.lat_first_deg + static_cast<double>(j0) * g.dlat_deg) * degree;
    const double lat1 = lat0 + g.dlat_deg * degree;
    const double x0 = std::sin(lat0), x1 = std::sin(lat1);
    const double x = std::sin(lat_deg * degree);
    const double t = x1 == x0 ? 0.0 : (x - x0) / (x1 - x0);
    return std::fma(t, f[index(g, column, j0 + 1)] - f[index(g, column, j0)],
                    f[index(g, column, j0)]);
}

struct State {
    double r2 = 0.0;
    double r1 = 0.0;
    double rho = 0.0;
};

struct AngularRectangle {
    double lon1 = 0.0;
    double lon2 = 0.0;
    double lat1 = 0.0;
    double lat2 = 0.0;
};

inline State row_state(const Grid& g, const std::vector<double>& r2,
                       const std::vector<double>& r1,
                       const std::vector<double>& rho, std::size_t row,
                       double lon_deg) {
    return {interpolate_row(g, r2, row, lon_deg),
            interpolate_row(g, r1, row, lon_deg),
            interpolate_row(g, rho, row, lon_deg)};
}

inline State column_state(const Grid& g, const std::vector<double>& r2,
                          const std::vector<double>& r1,
                          const std::vector<double>& rho, std::size_t column,
                          double lat_deg) {
    return {interpolate_column(g, r2, column, lat_deg),
            interpolate_column(g, r1, column, lat_deg),
            interpolate_column(g, rho, column, lat_deg)};
}

struct RowInterpolation {
    std::size_t i0 = 0;
    std::size_t i1 = 0;
    double t = 0.0;
};

inline RowInterpolation prepare_row_interpolation(const Grid& g,
                                                  double lon_deg) {
    double u = (lon_deg - g.lon_first_deg) / g.dlon_deg;
    if (g.global_lon) {
        u = std::fmod(u, static_cast<double>(g.nlon));
        if (u < 0.0)
            u += static_cast<double>(g.nlon);
        const std::size_t i0 = static_cast<std::size_t>(std::floor(u));
        return {i0, (i0 + 1) % g.nlon, u - static_cast<double>(i0)};
    }
    if (u <= 0.0)
        return {0, 0, 0.0};
    if (u >= static_cast<double>(g.nlon - 1))
        return {g.nlon - 1, g.nlon - 1, 0.0};
    const std::size_t i0 = static_cast<std::size_t>(std::floor(u));
    return {i0, i0 + 1, u - static_cast<double>(i0)};
}

struct ColumnInterpolation {
    std::size_t j0 = 0;
    std::size_t j1 = 0;
    double t = 0.0;
};

inline ColumnInterpolation prepare_column_interpolation(const Grid& g,
                                                        double lat_deg) {
    double u = (lat_deg - g.lat_first_deg) / g.dlat_deg;
    if (u <= 0.0)
        return {0, 0, 0.0};
    if (u >= static_cast<double>(g.nlat - 1))
        return {g.nlat - 1, g.nlat - 1, 0.0};
    const std::size_t j0 = static_cast<std::size_t>(std::floor(u));
    const double lat0 =
        (g.lat_first_deg + static_cast<double>(j0) * g.dlat_deg) * degree;
    const double lat1 = lat0 + g.dlat_deg * degree;
    const double x0 = std::sin(lat0), x1 = std::sin(lat1);
    const double x = std::sin(lat_deg * degree);
    return {j0, j0 + 1, x1 == x0 ? 0.0 : (x - x0) / (x1 - x0)};
}

inline double interpolate_prepared(const Grid& g,
                                   const std::vector<double>& field,
                                   const RowInterpolation& interpolation,
                                   std::size_t row) {
    const double a = field[index(g, interpolation.i0, row)];
    const double b = field[index(g, interpolation.i1, row)];
    return std::fma(interpolation.t, b - a, a);
}

inline double interpolate_prepared(const Grid& g,
                                   const std::vector<double>& field,
                                   const ColumnInterpolation& interpolation,
                                   std::size_t column) {
    const double a = field[index(g, column, interpolation.j0)];
    const double b = field[index(g, column, interpolation.j1)];
    return std::fma(interpolation.t, b - a, a);
}

inline State row_state(const Grid& g, const std::vector<double>& r2,
                       const std::vector<double>& r1,
                       const std::vector<double>& rho, std::size_t row,
                       const RowInterpolation& interpolation) {
    return {interpolate_prepared(g, r2, interpolation, row),
            interpolate_prepared(g, r1, interpolation, row),
            interpolate_prepared(g, rho, interpolation, row)};
}

inline State column_state(const Grid& g, const std::vector<double>& r2,
                          const std::vector<double>& r1,
                          const std::vector<double>& rho, std::size_t column,
                          const ColumnInterpolation& interpolation) {
    return {interpolate_prepared(g, r2, interpolation, column),
            interpolate_prepared(g, r1, interpolation, column),
            interpolate_prepared(g, rho, interpolation, column)};
}

inline std::size_t periodic_index(long long index_value,
                                  std::size_t count) {
    const long long n = static_cast<long long>(count);
    const long long wrapped = ((index_value % n) + n) % n;
    return static_cast<std::size_t>(wrapped);
}

inline double interpolate_row_smooth(const Grid& g,
                                     const std::vector<double>& field,
                                     std::size_t row, double lon_deg) {
    const double u = (lon_deg - g.lon_first_deg) / g.dlon_deg;
    const long long knot = static_cast<long long>(std::llround(u));
    const double half_width =
        std::min(0.05 / g.dlon_deg, 0.01);
    const bool valid = g.global_lon
        ? g.nlon >= 3
        : knot > 0 && knot + 1 < static_cast<long long>(g.nlon);
    if (!valid || std::abs(u - static_cast<double>(knot)) > half_width)
        return interpolate_row(g, field, row, lon_deg);

    auto value = [&](long long column) {
        const std::size_t i = g.global_lon
            ? periodic_index(column, g.nlon)
            : static_cast<std::size_t>(column);
        return field[index(g, i, row)];
    };
    const double fm = value(knot - 1);
    const double f0 = value(knot);
    const double fp = value(knot + 1);
    const double offset = u - static_cast<double>(knot);
    const double left = std::fma(offset, f0 - fm, f0);
    const double right = std::fma(offset, fp - f0, f0);
    const double weight = smoothstep5(
        (offset + half_width) / (2.0 * half_width));
    return std::fma(weight, right - left, left);
}

inline double interpolate_point_smooth(const Grid& g,
                                       const std::vector<double>& field,
                                       double lon_deg, double lat_deg) {
    const double u = (lat_deg - g.lat_first_deg) / g.dlat_deg;
    const long long knot = static_cast<long long>(std::llround(u));
    const double half_width =
        std::min(0.05 / g.dlat_deg, 0.01);
    if (knot > 0 && knot + 1 < static_cast<long long>(g.nlat) &&
        std::abs(u - static_cast<double>(knot)) <= half_width) {
        const std::size_t jm = static_cast<std::size_t>(knot - 1);
        const std::size_t j0 = static_cast<std::size_t>(knot);
        const std::size_t jp = static_cast<std::size_t>(knot + 1);
        const double fm = interpolate_row_smooth(g, field, jm, lon_deg);
        const double f0 = interpolate_row_smooth(g, field, j0, lon_deg);
        const double fp = interpolate_row_smooth(g, field, jp, lon_deg);
        const double xm = std::sin(
            (g.lat_first_deg + static_cast<double>(knot - 1) *
             g.dlat_deg) * degree);
        const double x0 = std::sin(
            (g.lat_first_deg + static_cast<double>(knot) *
             g.dlat_deg) * degree);
        const double xp = std::sin(
            (g.lat_first_deg + static_cast<double>(knot + 1) *
             g.dlat_deg) * degree);
        const double x = std::sin(lat_deg * degree);
        const double left = std::fma(
            (x - x0) / (x0 - xm), f0 - fm, f0);
        const double right = std::fma(
            (x - x0) / (xp - x0), fp - f0, f0);
        const double weight = smoothstep5(
            (u - static_cast<double>(knot) + half_width) /
            (2.0 * half_width));
        return std::fma(weight, right - left, left);
    }

    ColumnInterpolation lat = prepare_column_interpolation(g, lat_deg);
    const double south = interpolate_row_smooth(
        g, field, lat.j0, lon_deg);
    const double north = interpolate_row_smooth(
        g, field, lat.j1, lon_deg);
    return std::fma(lat.t, north - south, south);
}

inline State point_state_smooth(const Grid& g,
                                const std::vector<double>& r2,
                                const std::vector<double>& r1,
                                const std::vector<double>& rho,
                                double lon_deg, double lat_deg) {
    return {interpolate_point_smooth(g, r2, lon_deg, lat_deg),
            interpolate_point_smooth(g, r1, lon_deg, lat_deg),
            interpolate_point_smooth(g, rho, lon_deg, lat_deg)};
}


inline bool active(const State& s) {
    return std::isfinite(s.r1) && std::isfinite(s.r2) &&
           std::isfinite(s.rho) && s.r1 > 0.0 && s.r2 > s.r1 &&
           s.rho != 0.0;
}

template<class Prepared>
inline detail::Jet state_boundary_kernel(
    const State& state,const Prepared& prepared,
    const detail::DirectionCosine& cosine) {
    if (!active(state))
        return detail::Jet(0.0);
    Prism prism;
    prism.lat1 = -0.5 * PI;
    prism.lat2 = 0.5 * PI;
    prism.lon1 = -PI;
    prism.lon2 = PI;
    prism.r1 = state.r1;
    prism.r2 = state.r2;
    prism.density = 1.0;
    return (state.r2-state.r1)*detail::w_hat(prism,prepared.r,cosine);
}

template<class Prepared>
inline detail::Jet endpoint_u_difference(
    double positive_radius,double negative_radius,const Prepared& prepared,
    const detail::DirectionCosine& cosine) {
    if(positive_radius==negative_radius)
        return detail::Jet(0.0);
    const double lo=std::min(positive_radius,negative_radius);
    const double hi=std::max(positive_radius,negative_radius);
    Prism prism;
    prism.lat1=-0.5*PI; prism.lat2=0.5*PI;
    prism.lon1=-PI; prism.lon2=PI;
    prism.r1=lo; prism.r2=hi; prism.density=1.0;
    const detail::Jet difference=
        (hi-lo)*detail::w_hat(prism,prepared.r,cosine);
    return positive_radius>negative_radius?difference:-difference;
}

template<class Prepared>
inline detail::Jet material_jump_kernel(
    const State& positive,const State& negative,const Prepared& prepared,
    const detail::DirectionCosine& cosine) {
    if(!active(positive))
        return active(negative)
            ? -negative.rho*state_boundary_kernel(negative,prepared,cosine)
            : detail::Jet(0.0);
    if(!active(negative))
        return positive.rho*state_boundary_kernel(positive,prepared,cosine);

    const bool same_radii=
        positive.r1==negative.r1&&positive.r2==negative.r2;
    const bool same_density=positive.rho==negative.rho;
    if(same_radii&&same_density)
        return detail::Jet(0.0);
    if(same_radii)
        return (positive.rho-negative.rho)*
               state_boundary_kernel(positive,prepared,cosine);

    const detail::Jet delta_w=
        endpoint_u_difference(positive.r2,negative.r2,prepared,cosine)-
        endpoint_u_difference(positive.r1,negative.r1,prepared,cosine);
    if(same_density)
        return positive.rho*delta_w;

    const detail::Jet wp=state_boundary_kernel(positive,prepared,cosine);
    const detail::Jet wm=state_boundary_kernel(negative,prepared,cosine);
    const double rho_bar=0.5*(positive.rho+negative.rho);
    const double delta_rho=positive.rho-negative.rho;
    return rho_bar*delta_w+delta_rho*0.5*(wp+wm);
}

struct PreparedVectorGridLine {
    EdgeKind kind = EdgeKind::Latitude;
    double fixed = 0.0;
    detail::Jet r{};
    detail::JVec3 n{};
};

inline PreparedVectorGridLine prepare_vector_grid_line(
    EdgeKind kind, double fixed, const Point& point) {
    const detail::ScopedJetContext context(requested_component_context);
    const detail::JVec3 x = detail::observation_jet(point);
    const detail::Jet radius = detail::sqrt(detail::dot(x, x));
    return {kind, fixed, radius, detail::div(x, radius)};
}

inline Values prepared_jump_vector(const PreparedVectorGridLine& prepared,
                                   const State& positive,
                                   const State& negative,
                                   double parameter) {
    const detail::ScopedJetContext context(requested_component_context);
    Edge edge{prepared.kind, prepared.fixed, 0.0, 0.0, +1};
    detail::Vec3 source{}, tangent{};
    double alpha = 0.0;
    detail::source_geometry(edge, parameter, source, tangent, alpha);
    const detail::DirectionCosine cosine =
        detail::direction_cosine(prepared.n, source);
    const detail::Jet boundary =
        alpha * detail::dot(prepared.n, tangent) *
        material_jump_kernel(positive,negative,prepared,cosine);
    return detail::derivatives(boundary);
}

inline void accumulate_values(const Values& source, Values& destination) {
    for (std::size_t k = 0; k < destination.size(); ++k)
        destination[k] += source[k];
}

inline Values selected_latitude_aggregate_vector(
    const Grid& g, const std::vector<double>& r2,
    const std::vector<double>& r1, const std::vector<double>& rho,
    const std::vector<PreparedVectorGridLine>& prepared_lines,
    double longitude_rad, const std::vector<std::size_t>& selected) {
    const RowInterpolation interpolation =
        prepare_row_interpolation(g, longitude_rad / degree);
    Values total{}, correction{};
    for (const std::size_t j : selected) {
        const State south =
            j > 0 ? row_state(g, r2, r1, rho, j - 1, interpolation)
                  : State{};
        const State north =
            j < g.nlat ? row_state(g, r2, r1, rho, j, interpolation)
                       : State{};
        const Values value = prepared_jump_vector(
            prepared_lines[j], south, north, longitude_rad);
        for (std::size_t k = 0; k < total.size(); ++k)
            detail::compensated_add(value[k], total[k], correction[k]);
    }
    for (std::size_t k = 0; k < total.size(); ++k)
        total[k] += correction[k];
    return total;
}

inline Values selected_longitude_aggregate_vector(
    const Grid& g, const std::vector<double>& r2,
    const std::vector<double>& r1, const std::vector<double>& rho,
    const std::vector<PreparedVectorGridLine>& prepared_lines,
    double latitude_rad, const std::vector<std::size_t>& selected) {
    const ColumnInterpolation interpolation =
        prepare_column_interpolation(g, latitude_rad / degree);
    Values total{}, correction{};
    for (const std::size_t i : selected) {
        State west{}, east{};
        if (g.global_lon) {
            const std::size_t east_index = i % g.nlon;
            const std::size_t west_index =
                (east_index + g.nlon - 1) % g.nlon;
            west = column_state(g, r2, r1, rho, west_index,
                                interpolation);
            east = column_state(g, r2, r1, rho, east_index,
                                interpolation);
        } else {
            if (i > 0)
                west = column_state(g, r2, r1, rho, i - 1,
                                    interpolation);
            if (i < g.nlon)
                east = column_state(g, r2, r1, rho, i,
                                    interpolation);
        }
        const Values value = prepared_jump_vector(
            prepared_lines[i], west, east, latitude_rad);
        for (std::size_t k = 0; k < total.size(); ++k)
            detail::compensated_add(value[k], total[k], correction[k]);
    }
    for (std::size_t k = 0; k < total.size(); ++k)
        total[k] += correction[k];
    return total;
}

struct SharedPointCache {
    struct RadialCoefficients {
        double constant = 0.0;
        double radius_squared = 0.0;
        double inverse_radius = 0.0;
    };

    std::vector<PreparedVectorGridLine> latitude_lines;
    std::vector<PreparedVectorGridLine> longitude_lines;
    std::unordered_map<double, Values> latitude_all;
    std::unordered_map<double, Values> longitude_all;
    RadialCoefficients full_radial;
    bool has_full_radial = false;
};

inline std::vector<std::size_t> missing_indices(
    const std::vector<std::size_t>& all,
    const std::vector<std::size_t>& selected) {
    std::vector<std::size_t> missing;
    missing.reserve(all.size() - std::min(all.size(), selected.size()));
    std::size_t j = 0;
    for (const std::size_t value : all) {
        while (j < selected.size() && selected[j] < value)
            ++j;
        if (j == selected.size() || selected[j] != value)
            missing.push_back(value);
    }
    return missing;
}

using RadialCoefficients = SharedPointCache::RadialCoefficients;

inline void accumulate_radial_face(double radius, double weight,
                                   const Point& point, double sign,
                                   RadialCoefficients& coefficients) {
    const double signed_weight = sign * weight;
    if (point.r < radius) {
        coefficients.constant += signed_weight * radius * radius;
        coefficients.radius_squared -= signed_weight / 3.0;
    } else {
        coefficients.inverse_radius +=
            signed_weight * (2.0 / 3.0) * radius * radius * radius;
    }
}

inline Values radial_derivatives(const RadialCoefficients& coefficients,
                                const Point& point) {
    Values sum{};
    const detail::ScopedJetContext context(requested_component_context);
    const detail::JVec3 x = detail::observation_jet(point);
    const detail::Jet radius = detail::sqrt(detail::dot(x, x));
    const detail::Jet c0 = coefficients.constant +
        coefficients.radius_squared * radius * radius +
        coefficients.inverse_radius / radius;
    sum = detail::derivatives(c0);
    return sum;
}

template <class Function>
inline Values gauss1_vector(const Function& function, double a, double b) {
    Values value = function(0.5 * (a + b));
    for (double& item : value)
        item *= b - a;
    return value;
}

template <class Function>
inline Values gauss2_vector(const Function& function, double a, double b) {
    constexpr double node = 0.577350269189625764509148780501957;
    const double center = 0.5 * (a + b);
    const double half = 0.5 * (b - a);
    const Values left = function(center - half * node);
    const Values right = function(center + half * node);
    Values value{};
    for (std::size_t k = 0; k < value.size(); ++k)
        value[k] = half * (left[k] + right[k]);
    return value;
}

template <class Function>
inline Values gauss6_vector(const Function& function, double a, double b) {
    static constexpr std::array<double, 3> nodes{
        0.238619186083196908630501721681,
        0.661209386466264513661399595020,
        0.932469514203152027812301554494};
    static constexpr std::array<double, 3> weights{
        0.467913934572691047389870343990,
        0.360761573048138607569833513838,
        0.171324492379170345040296142173};
    const double center = 0.5 * (a + b);
    const double half = 0.5 * (b - a);
    Values sum{};
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double offset = half * nodes[i];
        const Values left = function(center - offset);
        const Values right = function(center + offset);
        for (std::size_t k = 0; k < sum.size(); ++k)
            sum[k] += weights[i] * (left[k] + right[k]);
    }
    for (double& value : sum)
        value *= half;
    return sum;
}

template <class Function>
inline Values gauss10_vector(const Function& function, double a, double b) {
    static constexpr std::array<double, 5> nodes{
        0.148874338981631210884826001130,
        0.433395394129247190799265943166,
        0.679409568299024406234327365115,
        0.865063366688984510732096688423,
        0.973906528517171720077964012084};
    static constexpr std::array<double, 5> weights{
        0.295524224714752870173892994651,
        0.269266719309996355091226921569,
        0.219086362515982043995534934228,
        0.149451349150580593145776339658,
        0.066671344308688137593568809893};
    const double center = 0.5 * (a + b);
    const double half = 0.5 * (b - a);
    Values sum{};
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double offset = half * nodes[i];
        const Values left = function(center - offset);
        const Values right = function(center + offset);
        for (std::size_t k = 0; k < sum.size(); ++k)
            sum[k] += weights[i] * (left[k] + right[k]);
    }
    for (double& value : sum)
        value *= half;
    return sum;
}



template <class Function>
inline Values integrate_vector_line(
    const Function& function, const Edge& edge, const Point& point,
    double interpolation_origin, double interpolation_step) {
    if (!(interpolation_step > 0.0))
        throw std::invalid_argument(
            "TLI-MMSI panel integration requires a positive grid step");

    const double geometric = tli_mmsi_detail::clamp_value(
        detail::geometric_projection(edge, point), edge.a, edge.b);

    std::vector<double> panel{edge.a, geometric, edge.b};
    const long long first = static_cast<long long>(std::ceil(
        (edge.a - interpolation_origin) / interpolation_step));
    const long long last = static_cast<long long>(std::floor(
        (edge.b - interpolation_origin) / interpolation_step));
    for (long long i = first; i <= last; ++i) {
        const double knot = interpolation_origin +
                            static_cast<double>(i) * interpolation_step;
        if (knot > edge.a && knot < edge.b)
            panel.push_back(knot);
    }
    std::sort(panel.begin(), panel.end());
    const double position_tolerance =
        128.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(edge.b - edge.a));
    panel.erase(std::unique(
        panel.begin(), panel.end(), [&](double a, double b) {
            return std::abs(a - b) <= position_tolerance;
        }), panel.end());

    const bool fine_grid = interpolation_step <= 0.125 * degree;
    constexpr double accurate_window = 5.0 * degree;
    auto projection_at = [&](double coordinate) {
        return std::abs(coordinate - geometric) <= position_tolerance;
    };

    Values total{}, correction{};
    auto accumulate = [&](const Values& value) {
        for (std::size_t k = 0; k < total.size(); ++k)
            detail::compensated_add(value[k], total[k], correction[k]);
    };
    for (std::size_t i = 1; i < panel.size(); ++i) {
        const double a = panel[i - 1], b = panel[i];
        const bool projection_left = projection_at(a);
        const bool projection_right = projection_at(b);
        if (projection_left || projection_right) {
            auto mapped = [&](double s) {
                double coordinate = 0.0, jacobian = 0.0;
                if (projection_left && projection_right) {
                    coordinate = a + 0.5 * (b - a) *
                                         (1.0 - std::cos(s));
                    jacobian = 0.5 * (b - a) * std::sin(s);
                } else if (projection_left) {
                    coordinate = a + (b - a) * (1.0 - std::cos(s));
                    jacobian = (b - a) * std::sin(s);
                } else {
                    coordinate = a + (b - a) * std::sin(s);
                    jacobian = (b - a) * std::cos(s);
                }
                Values value = function(coordinate);
                for (double& item : value)
                    item *= jacobian;
                return value;
            };
            accumulate(gauss10_vector(
                mapped, 0.0,
                projection_left && projection_right ? PI : 0.5 * PI));
            continue;
        }

        if (!fine_grid) {
            accumulate(gauss10_vector(function, a, b));
            continue;
        }
        double distance = std::abs(0.5 * (a + b) - geometric);
        if (edge.kind == EdgeKind::Latitude)
            distance = std::abs(std::remainder(distance, 2.0 * PI));
        const bool accurate = distance <= accurate_window;
        accumulate(accurate ? gauss2_vector(function, a, b)
                            : gauss1_vector(function, a, b));
    }
    for (std::size_t k = 0; k < total.size(); ++k)
        total[k] += correction[k];
    return total;
}

inline void add_scaled(const Values& source, double scale, Values& sum) {
    for (std::size_t k = 0; k < sum.size(); ++k)
        sum[k] = std::fma(scale, source[k], sum[k]);
}

inline Values explicit_sum_outside_rectangle(
    const Grid& g, const std::vector<double>& r2,
    const std::vector<double>& r1, const std::vector<double>& rho,
    const Point& point, const AngularRectangle& hole) {
    double constant_coefficient = 0.0;
    double radius_squared_coefficient = 0.0;
    double inverse_radius_coefficient = 0.0;
    const double dlon = g.dlon_deg * degree;
    const double lon_west = g.lon_west_deg() * degree;
    const double lat_south = g.lat_south_deg() * degree;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : constant_coefficient, radius_squared_coefficient, inverse_radius_coefficient) schedule(static) if(!omp_in_parallel())
#endif
    for (long long jj = 0; jj < static_cast<long long>(g.nlat); ++jj) {
        const std::size_t j = static_cast<std::size_t>(jj);
        const double lat1 = std::max(-0.5 * PI,
            lat_south + static_cast<double>(j) * g.dlat_deg * degree);
        const double lat2 = std::min(0.5 * PI,
            lat_south + static_cast<double>(j + 1) * g.dlat_deg * degree);
        for (std::size_t i = 0; i < g.nlon; ++i) {
            const std::size_t p = index(g, i, j);
            const State state{r2[p], r1[p], rho[p]};
            if (!active(state))
                continue;
            const double lon1 = lon_west + static_cast<double>(i) * dlon;
            const double lon2 = lon1 + dlon;
            const double overlap_lon = std::max(
                0.0, std::min(lon2, hole.lon2) - std::max(lon1, hole.lon1));
            const double overlap_lat1 = std::max(lat1, hole.lat1);
            const double overlap_lat2 = std::min(lat2, hole.lat2);
            const double overlap = overlap_lon > 0.0 && overlap_lat2 > overlap_lat1
                ? overlap_lon * (std::sin(overlap_lat2) - std::sin(overlap_lat1))
                : 0.0;
            const double full = dlon * (std::sin(lat2) - std::sin(lat1));
            const double weight = std::max(0.0, full - overlap) * state.rho;
            auto accumulate_face = [&](double a, double sign) {
                const double signed_weight = sign * weight;
                if (point.r < a) {
                    constant_coefficient += signed_weight * a * a;
                    radius_squared_coefficient -= signed_weight / 3.0;
                } else {
                    inverse_radius_coefficient +=
                        signed_weight * (2.0 / 3.0) * a * a * a;
                }
            };
            accumulate_face(state.r2, +1.0);
            accumulate_face(state.r1, -1.0);
        }
    }
    Values sum{};
    {
        const detail::ScopedJetContext context(requested_component_context);
        const detail::JVec3 x = detail::observation_jet(point);
        const detail::Jet radius = detail::sqrt(detail::dot(x, x));
        const detail::Jet c0 = constant_coefficient +
            radius_squared_coefficient * radius * radius +
            inverse_radius_coefficient / radius;
        sum = detail::derivatives(c0);
    }
    return sum;
}

inline Values explicit_sum_outside_rectangle_shared(
    const Grid& g, const std::vector<double>& r2,
    const std::vector<double>& r1, const std::vector<double>& rho,
    const Point& point, const AngularRectangle& hole,
    SharedPointCache& cache) {
    const double dlon = g.dlon_deg * degree;
    const double dlat = g.dlat_deg * degree;
    const double lon_west = g.lon_west_deg() * degree;
    const double lat_south = g.lat_south_deg() * degree;
    if (!cache.has_full_radial) {
        const double full_weight_scale = dlon;
        double full_constant = 0.0;
        double full_radius_squared = 0.0;
        double full_inverse_radius = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : full_constant, full_radius_squared, full_inverse_radius) schedule(static) if(!omp_in_parallel())
#endif
        for (long long jj = 0; jj < static_cast<long long>(g.nlat); ++jj) {
            const std::size_t j = static_cast<std::size_t>(jj);
            const double lat1 = std::max(-0.5 * PI,
                lat_south + static_cast<double>(j) * dlat);
            const double lat2 = std::min(0.5 * PI,
                lat_south + static_cast<double>(j + 1) * dlat);
            for (std::size_t i = 0; i < g.nlon; ++i) {
                const std::size_t p = index(g, i, j);
                const State state{r2[p], r1[p], rho[p]};
                if (!active(state))
                    continue;
                const double weight = full_weight_scale *
                    (std::sin(lat2) - std::sin(lat1)) * state.rho;
                RadialCoefficients local{};
                accumulate_radial_face(state.r2, weight, point, +1.0, local);
                accumulate_radial_face(state.r1, weight, point, -1.0, local);
                full_constant += local.constant;
                full_radius_squared += local.radius_squared;
                full_inverse_radius += local.inverse_radius;
            }
        }
        cache.full_radial = {full_constant, full_radius_squared,
                             full_inverse_radius};
        cache.has_full_radial = true;
    }

    RadialCoefficients coefficients = cache.full_radial;
    const long long first_lat = std::max<long long>(
        0, static_cast<long long>(std::floor((hole.lat1 - lat_south) / dlat)) - 1);
    const long long last_lat = std::min<long long>(
        static_cast<long long>(g.nlat) - 1,
        static_cast<long long>(std::ceil((hole.lat2 - lat_south) / dlat)) + 1);
    const long long first_lon = std::max<long long>(
        0, static_cast<long long>(std::floor((hole.lon1 - lon_west) / dlon)) - 1);
    const long long last_lon = std::min<long long>(
        static_cast<long long>(g.nlon) - 1,
        static_cast<long long>(std::ceil((hole.lon2 - lon_west) / dlon)) + 1);
    for (long long jj = first_lat; jj <= last_lat; ++jj) {
        const std::size_t j = static_cast<std::size_t>(jj);
        const double lat1 = std::max(-0.5 * PI,
            lat_south + static_cast<double>(j) * dlat);
        const double lat2 = std::min(0.5 * PI,
            lat_south + static_cast<double>(j + 1) * dlat);
        const double overlap_lat1 = std::max(lat1, hole.lat1);
        const double overlap_lat2 = std::min(lat2, hole.lat2);
        if (!(overlap_lat2 > overlap_lat1))
            continue;
        for (long long ii = first_lon; ii <= last_lon; ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            const std::size_t p = index(g, i, j);
            const State state{r2[p], r1[p], rho[p]};
            if (!active(state))
                continue;
            const double lon1 = lon_west + static_cast<double>(i) * dlon;
            const double lon2 = lon1 + dlon;
            const double overlap_lon = std::max(
                0.0, std::min(lon2, hole.lon2) -
                     std::max(lon1, hole.lon1));
            if (!(overlap_lon > 0.0))
                continue;
            const double overlap = overlap_lon *
                (std::sin(overlap_lat2) - std::sin(overlap_lat1)) * state.rho;
            accumulate_radial_face(state.r2, overlap, point, -1.0,
                                   coefficients);
            accumulate_radial_face(state.r1, overlap, point, +1.0,
                                   coefficients);
        }
    }
    return radial_derivatives(coefficients, point);
}

inline bool coordinate_inside(double value, double lower, double upper) {
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(value), std::abs(lower), std::abs(upper)});
    return value >= lower - tolerance && value <= upper + tolerance;
}

inline long long cell_below(double coordinate, double lower, double step,
                            std::size_t count) {
    if (coordinate <= lower)
        return -1;
    double u = (coordinate - lower) / step;
    const double nearest = std::round(u);
    if (std::abs(u - nearest) <=
        128.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(u)))
        u = nearest;
    return std::min<long long>(static_cast<long long>(count) - 1,
                               static_cast<long long>(std::ceil(u)) - 1);
}

inline long long cell_above(double coordinate, double lower, double step,
                            std::size_t count) {
    const double upper = lower + static_cast<double>(count) * step;
    if (coordinate >= upper)
        return -1;
    double u = (coordinate - lower) / step;
    const double nearest = std::round(u);
    if (std::abs(u - nearest) <=
        128.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(u)))
        u = nearest;
    return std::max<long long>(0,
        static_cast<long long>(std::floor(u)));
}

template <class Function>
inline void add_gauss6_interval(const Function& function, double a, double b,
                                Values& total) {
    if (b > a)
        add_scaled(gauss6_vector(function, a, b), 1.0, total);
}

template <class Function>
inline void add_composite_gauss6(const Function& function, double a, double b,
                                 double maximum_width, Values& total) {
    if (!(b > a))
        return;
    const std::size_t pieces = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil((b - a) / maximum_width)));
    for (std::size_t i = 0; i < pieces; ++i) {
        const double left = std::fma(
            static_cast<double>(i) / static_cast<double>(pieces), b - a, a);
        const double right = i + 1 == pieces ? b : std::fma(
            static_cast<double>(i + 1) / static_cast<double>(pieces),
            b - a, a);
        add_gauss6_interval(function, left, right, total);
    }
}

template <class Function>
inline void add_accurate_interval(const Function& function, EdgeKind kind,
                                  double fixed, double a, double b,
                                  const Point& point,
                                  double interpolation_origin,
                                  double interpolation_step, Values& total) {
    if (!(b > a))
        return;
    const Edge edge{kind, fixed, a, b, +1};
    add_scaled(integrate_vector_line(function, edge, point,
                                     interpolation_origin,
                                     interpolation_step),
               1.0, total);
}

template <class Function>
inline void add_accuracy_correction(const Function& function, EdgeKind kind,
                                    double fixed, double a, double b,
                                    const Point& point,
                                    double interpolation_origin,
                                    double interpolation_step, double weight,
                                    Values& total) {
    if (!(b > a) || weight == 0.0)
        return;
    const Edge edge{kind, fixed, a, b, +1};
    const Values accurate = integrate_vector_line(
        function, edge, point, interpolation_origin, interpolation_step);
    const Values baseline = gauss6_vector(function, a, b);
    add_scaled(accurate, weight, total);
    add_scaled(baseline, -weight, total);
}

template <class Function>
inline void add_line_outside_hole(
    const Function& function, EdgeKind kind, double fixed,
    double domain_a, double domain_b, double local_a, double local_b,
    double hole_a, double hole_b, const Point& point,
    double interpolation_origin, double interpolation_step,
    double far_maximum_width, Values& total) {
    auto add_piece = [&](double a, double b) {
        if (!(b > a))
            return;
        add_composite_gauss6(function, a, std::min(b, local_a),
                             far_maximum_width, total);
        add_accurate_interval(function, kind, fixed,
                              std::max(a, local_a), std::min(b, local_b),
                              point, interpolation_origin,
                              interpolation_step, total);
        add_composite_gauss6(function, std::max(a, local_b), b,
                             far_maximum_width, total);
    };
    add_piece(domain_a, std::min(domain_b, hole_a));
    add_piece(std::max(domain_a, hole_b), domain_b);
}

inline Values constant_tesseroid_estimate(
    const AngularRectangle& rectangle, const State& state, const Point& point,
    double absolute_tolerance, double relative_tolerance) {
    Values values{};
    if (!active(state) || rectangle.lon2 <= rectangle.lon1 ||
        rectangle.lat2 <= rectangle.lat1)
        return values;
    Prism prism;
    prism.lon1 = rectangle.lon1;
    prism.lon2 = rectangle.lon2;
    prism.lat1 = rectangle.lat1;
    prism.lat2 = rectangle.lat2;
    prism.r1 = state.r1;
    prism.r2 = state.r2;
    prism.density = state.rho;
    EstimateOptions options;
    options.gravitational_constant = 1.0;
    constexpr double tolerance_floor = 1.0e-10;
    options.integration.atol = std::max(
        std::max(absolute_tolerance, 0.0), tolerance_floor);
    options.integration.rtol = std::max(
        std::max(relative_tolerance, 0.0), tolerance_floor);
    for (std::size_t k = 0; k < values.size(); ++k) {
        if ((requested_component_mask & (1u << k)) == 0)
            continue;
        values[k] = estimate(static_cast<Component>(k), prism, point, options);
    }
    return values;
}

inline Values point_estimate_configuration(
    const Grid& g, const std::vector<double>& r2,
    const std::vector<double>& r1, const std::vector<double>& rho,
    double observation_lat_deg, double observation_lon_deg,
    double observation_radius, double absolute_tolerance,
    double relative_tolerance, double hole_lat_center_deg,
    double hole_lon_center_deg, std::size_t component_context,
    std::uint32_t canonical_component_mask,
    SharedPointCache* shared_cache = nullptr) {
    const ScopedComponentContext requested_context(
        component_context, canonical_component_mask);
    const Point point{observation_lat_deg * degree,
                      detail::normalize_longitude(observation_lon_deg *
                                                  degree),
                      observation_radius};
    const double observation_lon_normalized_deg = point.lon / degree;
    const double lon_a = g.lon_west_deg() * degree;
    const double lon_b = g.lon_east_deg() * degree;
    double lon_projection = point.lon;
    while (lon_projection < lon_a)
        lon_projection += 2.0 * PI;
    while (lon_projection > lon_b)
        lon_projection -= 2.0 * PI;
    double hole_lon_projection =
        detail::normalize_longitude(hole_lon_center_deg * degree);
    while (hole_lon_projection < lon_a)
        hole_lon_projection += 2.0 * PI;
    while (hole_lon_projection > lon_b)
        hole_lon_projection -= 2.0 * PI;
    const double lat_a = g.lat_south_deg() * degree;
    const double lat_b = g.lat_north_deg() * degree;
    const double hole_lat_center = hole_lat_center_deg * degree;

    const AngularRectangle hole{
        std::max(lon_a, hole_lon_projection - 0.5 * g.dlon_deg * degree),
        std::min(lon_b, hole_lon_projection + 0.5 * g.dlon_deg * degree),
        std::max(lat_a, hole_lat_center - 0.5 * g.dlat_deg * degree),
        std::min(lat_b, hole_lat_center + 0.5 * g.dlat_deg * degree)};
    const bool fine_grid =
        std::max(g.dlon_deg, g.dlat_deg) < 1.0;
    const double local_lon_half_width = degree;
    const double local_lat_half_width = degree;
    const AngularRectangle local{
        std::max(lon_a, lon_projection - local_lon_half_width),
        std::min(lon_b, lon_projection + local_lon_half_width),
        std::max(lat_a, point.lat - local_lat_half_width),
        std::min(lat_b, point.lat + local_lat_half_width)};

    Values total = shared_cache
        ? explicit_sum_outside_rectangle_shared(
            g, r2, r1, rho, point, hole, *shared_cache)
        : explicit_sum_outside_rectangle(g, r2, r1, rho, point, hole);
    const std::size_t latitude_first =
        g.lat_south_deg() == -90.0 ? 1U : 0U;
    const std::size_t latitude_last =
        g.lat_north_deg() == 90.0 ? g.nlat - 1 : g.nlat;
    std::vector<PreparedVectorGridLine> local_latitude_lines;
    if (shared_cache) {
        if (shared_cache->latitude_lines.empty()) {
            shared_cache->latitude_lines.resize(g.nlat + 1);
            for (std::size_t j = 0; j <= g.nlat; ++j) {
                const double fixed = lat_a + static_cast<double>(j) *
                                               g.dlat_deg * degree;
                shared_cache->latitude_lines[j] =
                    prepare_vector_grid_line(EdgeKind::Latitude, fixed, point);
            }
        }
    } else {
        local_latitude_lines.resize(g.nlat + 1);
        for (std::size_t j = 0; j <= g.nlat; ++j) {
            const double fixed = lat_a + static_cast<double>(j) *
                                           g.dlat_deg * degree;
            local_latitude_lines[j] =
                prepare_vector_grid_line(EdgeKind::Latitude, fixed, point);
        }
    }
    const auto& latitude_lines = shared_cache
        ? shared_cache->latitude_lines : local_latitude_lines;
    std::vector<std::size_t> all_latitude;
    std::vector<std::size_t> far_latitude;
    if (latitude_last >= latitude_first)
        for (std::size_t j = latitude_first; j <= latitude_last; ++j) {
            all_latitude.push_back(j);
            const double fixed = latitude_lines[j].fixed;
            if (!coordinate_inside(fixed, hole.lat1, hole.lat2))
                far_latitude.push_back(j);
        }
    auto latitude_function = [&](const std::vector<std::size_t>& selected) {
        const auto* selected_lines = &selected;
        return [&, selected_lines](double longitude) {
            if (!shared_cache || selected_lines->size() <= 1)
                return selected_latitude_aggregate_vector(
                    g, r2, r1, rho, latitude_lines, longitude,
                    *selected_lines);
            const auto cached = shared_cache->latitude_all.find(longitude);
            Values total{};
            if (cached == shared_cache->latitude_all.end()) {
                total = selected_latitude_aggregate_vector(
                    g, r2, r1, rho, latitude_lines, longitude,
                    all_latitude);
                shared_cache->latitude_all.emplace(longitude, total);
            } else {
                total = cached->second;
            }
            const std::vector<std::size_t> omitted = missing_indices(
                all_latitude, *selected_lines);
            if (!omitted.empty()) {
                const Values removed = selected_latitude_aggregate_vector(
                    g, r2, r1, rho, latitude_lines, longitude, omitted);
                for (std::size_t k = 0; k < total.size(); ++k)
                    total[k] -= removed[k];
            }
            return total;
        };
    };
    const double far_width = fine_grid ? 60.0 * degree : 10.0 * degree;
    const double switching_far_width = fine_grid ? far_width : degree;
    if (!far_latitude.empty()) {
        const auto function = latitude_function(far_latitude);
        add_composite_gauss6(function, lon_a, lon_b,
                             far_width, total);
    }
    for (const std::size_t j : all_latitude) {
        const std::vector<std::size_t> one{j};
        const auto single = latitude_function(one);
        const double fixed = latitude_lines[j].fixed;
        const bool crosses_hole = coordinate_inside(fixed, hole.lat1,
                                                     hole.lat2);
        if (crosses_hole) {
            add_line_outside_hole(
                single, EdgeKind::Latitude, fixed, lon_a, lon_b,
                local.lon1, local.lon2, hole.lon1, hole.lon2, point,
                g.lon_first_deg * degree, g.dlon_deg * degree,
                switching_far_width, total);
        }
        if (!fine_grid && !crosses_hole &&
            std::min(std::abs(fixed - hole.lat1),
                     std::abs(fixed - hole.lat2)) <= g.dlat_deg * degree) {
            Values fine{}, coarse{};
            add_composite_gauss6(single, lon_a, lon_b,
                                 switching_far_width, fine);
            add_composite_gauss6(single, lon_a, lon_b, far_width, coarse);
            accumulate_values(fine, total);
            add_scaled(coarse, -1.0, total);
        }
        if (fine_grid && !crosses_hole) {
            const double distance = std::abs(fixed - point.lat);
            const double transition = g.dlat_deg * degree;
            const double weight = clamp_value(
                (degree + transition - distance) / transition, 0.0, 1.0);
            add_accuracy_correction(
                single, EdgeKind::Latitude, fixed, local.lon1,
                local.lon2, point, g.lon_first_deg * degree,
                g.dlon_deg * degree, weight, total);
        }
    }

    const std::size_t longitude_count =
        g.global_lon ? g.nlon : g.nlon + 1;
    std::vector<PreparedVectorGridLine> local_longitude_lines;
    if (shared_cache) {
        if (shared_cache->longitude_lines.empty()) {
            shared_cache->longitude_lines.resize(longitude_count);
            for (std::size_t i = 0; i < longitude_count; ++i) {
                const double fixed = lon_a + static_cast<double>(i) *
                                               g.dlon_deg * degree;
                shared_cache->longitude_lines[i] =
                    prepare_vector_grid_line(EdgeKind::Longitude, fixed, point);
            }
        }
    } else {
        local_longitude_lines.resize(longitude_count);
        for (std::size_t i = 0; i < longitude_count; ++i) {
            const double fixed = lon_a + static_cast<double>(i) *
                                           g.dlon_deg * degree;
            local_longitude_lines[i] =
                prepare_vector_grid_line(EdgeKind::Longitude, fixed, point);
        }
    }
    const auto& longitude_lines = shared_cache
        ? shared_cache->longitude_lines : local_longitude_lines;
    std::vector<std::size_t> all_longitude(longitude_count);
    std::vector<std::size_t> far_longitude;
    for (std::size_t i = 0; i < longitude_count; ++i) {
        all_longitude[i] = i;
        if (!coordinate_inside(longitude_lines[i].fixed,
                               hole.lon1, hole.lon2))
            far_longitude.push_back(i);
    }
    auto longitude_function = [&](const std::vector<std::size_t>& selected) {
        const auto* selected_lines = &selected;
        return [&, selected_lines](double latitude) {
            if (!shared_cache || selected_lines->size() <= 1)
                return selected_longitude_aggregate_vector(
                    g, r2, r1, rho, longitude_lines, latitude,
                    *selected_lines);
            const auto cached = shared_cache->longitude_all.find(latitude);
            Values total{};
            if (cached == shared_cache->longitude_all.end()) {
                total = selected_longitude_aggregate_vector(
                    g, r2, r1, rho, longitude_lines, latitude,
                    all_longitude);
                shared_cache->longitude_all.emplace(latitude, total);
            } else {
                total = cached->second;
            }
            const std::vector<std::size_t> omitted = missing_indices(
                all_longitude, *selected_lines);
            if (!omitted.empty()) {
                const Values removed = selected_longitude_aggregate_vector(
                    g, r2, r1, rho, longitude_lines, latitude, omitted);
                for (std::size_t k = 0; k < total.size(); ++k)
                    total[k] -= removed[k];
            }
            return total;
        };
    };
    if (!far_longitude.empty()) {
        const auto function = longitude_function(far_longitude);
        add_composite_gauss6(function, lat_a, lat_b,
                             far_width, total);
    }
    for (const std::size_t i : all_longitude) {
        const std::vector<std::size_t> one{i};
        const auto single = longitude_function(one);
        const double fixed = longitude_lines[i].fixed;
        const bool crosses_hole = coordinate_inside(fixed, hole.lon1,
                                                     hole.lon2);
        if (crosses_hole) {
            add_line_outside_hole(
                single, EdgeKind::Longitude, fixed, lat_a, lat_b,
                local.lat1, local.lat2, hole.lat1, hole.lat2, point,
                g.lat_first_deg * degree, g.dlat_deg * degree,
                switching_far_width, total);
        }
        if (!fine_grid && !crosses_hole &&
            std::min(std::abs(fixed - hole.lon1),
                     std::abs(fixed - hole.lon2)) <= g.dlon_deg * degree) {
            Values fine{}, coarse{};
            add_composite_gauss6(single, lat_a, lat_b,
                                 switching_far_width, fine);
            add_composite_gauss6(single, lat_a, lat_b, far_width, coarse);
            accumulate_values(fine, total);
            add_scaled(coarse, -1.0, total);
        }
        if (fine_grid && !crosses_hole) {
            const double distance = g.global_lon
                ? std::abs(std::remainder(fixed - point.lon, 2.0 * PI))
                : std::abs(fixed - point.lon);
            const double transition = g.dlon_deg * degree;
            const double weight = clamp_value(
                (degree + transition - distance) / transition, 0.0, 1.0);
            add_accuracy_correction(
                single, EdgeKind::Longitude, fixed, local.lat1,
                local.lat2, point, g.lat_first_deg * degree,
                g.dlat_deg * degree, weight, total);
        }
    }

    const double lat_step = g.dlat_deg * degree;
    const double lon_step = g.dlon_deg * degree;
    const long long south_row = cell_below(hole.lat1, lat_a, lat_step,
                                           g.nlat);
    const long long north_row = cell_above(hole.lat2, lat_a, lat_step,
                                           g.nlat);
    auto add_hole_latitude_edge = [&](double fixed, long long row,
                                      bool outside_positive) {
        const PreparedVectorGridLine prepared =
            prepare_vector_grid_line(EdgeKind::Latitude, fixed, point);
        const auto function = [&](double longitude) {
            const State outside = row >= 0
                ? row_state(g, r2, r1, rho, static_cast<std::size_t>(row),
                            longitude / degree)
                : State{};
            return outside_positive
                ? prepared_jump_vector(prepared, outside, State{}, longitude)
                : prepared_jump_vector(prepared, State{}, outside, longitude);
        };
        add_accurate_interval(function, EdgeKind::Latitude, fixed,
                              hole.lon1, hole.lon2, point,
                              g.lon_first_deg * degree, lon_step, total);
    };
    add_hole_latitude_edge(hole.lat1, south_row, true);
    add_hole_latitude_edge(hole.lat2, north_row, false);

    const long long west_column = cell_below(hole.lon1, lon_a, lon_step,
                                             g.nlon);
    const long long east_column = cell_above(hole.lon2, lon_a, lon_step,
                                             g.nlon);
    auto add_hole_longitude_edge = [&](double fixed, long long column,
                                       bool outside_positive) {
        const PreparedVectorGridLine prepared =
            prepare_vector_grid_line(EdgeKind::Longitude, fixed, point);
        const auto function = [&](double latitude) {
            const State outside = column >= 0
                ? column_state(g, r2, r1, rho,
                               static_cast<std::size_t>(column),
                               latitude / degree)
                : State{};
            return outside_positive
                ? prepared_jump_vector(prepared, outside, State{}, latitude)
                : prepared_jump_vector(prepared, State{}, outside, latitude);
        };
        add_accurate_interval(function, EdgeKind::Longitude, fixed,
                              hole.lat1, hole.lat2, point,
                              g.lat_first_deg * degree, lat_step, total);
    };
    add_hole_longitude_edge(hole.lon1, west_column, true);
    add_hole_longitude_edge(hole.lon2, east_column, false);

    for (double& value : total)
        value *= 0.5;
    const State local_state = point_state_smooth(
        g, r2, r1, rho, observation_lon_normalized_deg,
        observation_lat_deg);
    const Values local_tesseroid = constant_tesseroid_estimate(
        hole, local_state, point, absolute_tolerance, relative_tolerance);
    accumulate_values(local_tesseroid, total);
    return total;
}

struct BlendAxis {
    std::array<double, 2> center{};
    std::array<double, 2> weight{1.0, 0.0};
    std::size_t count = 1;
};

inline BlendAxis blend_axis(double coordinate, double first, double step,
                            std::size_t count, bool periodic) {
    BlendAxis axis;
    axis.center[0] = coordinate;
    if (!(step > 0.0) || count == 0)
        return axis;

    const double u = (coordinate - first) / step;
    const long long nearest = static_cast<long long>(std::llround(u));
    if (!periodic &&
        (nearest < 0 || nearest >= static_cast<long long>(count)))
        return axis;

    const double knot = first + static_cast<double>(nearest) * step;
    const double distance = periodic
        ? std::remainder(coordinate - knot, 360.0)
        : coordinate - knot;
    const double half_width = std::min(0.05, 0.01 * step);
    if (!(half_width > 0.0) || std::abs(distance) > half_width)
        return axis;

    const double w = smoothstep5(
        (distance + half_width) / (2.0 * half_width));
    axis.center = {knot - half_width, knot + half_width};
    axis.weight = {1.0 - w, w};
    axis.count = 2;
    return axis;
}

inline Values point_estimate_tli(
    const Grid& g, const std::vector<double>& r2,
    const std::vector<double>& r1, const std::vector<double>& rho,
    double observation_lat_deg, double observation_lon_deg,
    double observation_radius, double absolute_tolerance,
    double relative_tolerance, std::size_t component_context = 20,
    std::uint32_t canonical_component_mask = (1u << 20) - 1u) {
    const double longitude = g.global_lon
        ? detail::normalize_longitude(observation_lon_deg * degree) / degree
        : observation_lon_deg;
    const BlendAxis longitude_axis = blend_axis(
        longitude, g.lon_first_deg, g.dlon_deg, g.nlon, g.global_lon);
    const BlendAxis latitude_axis = blend_axis(
        observation_lat_deg, g.lat_first_deg, g.dlat_deg, g.nlat, false);

    SharedPointCache shared_cache;
    SharedPointCache* shared_cache_pointer =
        latitude_axis.count * longitude_axis.count > 1
            ? &shared_cache : nullptr;
    Values total{};
    for (std::size_t j = 0; j < latitude_axis.count; ++j)
        for (std::size_t i = 0; i < longitude_axis.count; ++i) {
            const double weight =
                latitude_axis.weight[j] * longitude_axis.weight[i];
            const Values value = point_estimate_configuration(
                g, r2, r1, rho, observation_lat_deg, observation_lon_deg,
                observation_radius, absolute_tolerance, relative_tolerance,
                latitude_axis.center[j], longitude_axis.center[i],
                component_context, canonical_component_mask,
                shared_cache_pointer);
            add_scaled(value, weight, total);
        }
    return total;
}

inline State point_state(const Grid& g, const std::vector<double>& r2,
                         const std::vector<double>& r1,
                         const std::vector<double>& rho,
                         double lon_deg, double lat_deg) {
    const ColumnInterpolation latitude =
        prepare_column_interpolation(g, lat_deg);
    auto interpolate = [&](const std::vector<double>& field) {
        const double south = interpolate_row(
            g, field, latitude.j0, lon_deg);
        const double north = interpolate_row(
            g, field, latitude.j1, lon_deg);
        return std::fma(latitude.t, north - south, south);
    };
    return {interpolate(r2), interpolate(r1), interpolate(rho)};
}

struct HybridNearModel {
    Grid grid;
    std::vector<double> r2;
    std::vector<double> r1;
    std::vector<double> rho;
    AngularRectangle rectangle;
};

inline HybridNearModel make_hybrid_near_model(
    const Grid& source_grid, const std::vector<double>& source_r2,
    const std::vector<double>& source_r1,
    const std::vector<double>& source_rho, double observation_lon_deg,
    double observation_lat_deg) {
    const double lon_half_width_deg =
        std::max(1.0, 0.5 * source_grid.dlon_deg);
    const double lat_half_width_deg =
        std::max(1.0, 0.5 * source_grid.dlat_deg);
    const std::size_t nlon = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::llround(
               2.0 * lon_half_width_deg / source_grid.dlon_deg)));
    const double lon_first = observation_lon_deg - lon_half_width_deg +
                             0.5 * source_grid.dlon_deg;
    const double lat_south = std::max(
        -90.0, observation_lat_deg - lat_half_width_deg);
    const double lat_north = std::min(
        90.0, observation_lat_deg + lat_half_width_deg);
    const double lat_span = lat_north - lat_south;
    const std::size_t clipped_nlat = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::llround(
               lat_span / source_grid.dlat_deg)));
    const double lat_first = lat_south + 0.5 * source_grid.dlat_deg;

    HybridNearModel model;
    model.grid.lon_first_deg = lon_first;
    model.grid.lon_last_deg = lon_first +
        static_cast<double>(nlon - 1) * source_grid.dlon_deg;
    model.grid.lat_first_deg = lat_first;
    model.grid.lat_last_deg = lat_first +
        static_cast<double>(clipped_nlat - 1) * source_grid.dlat_deg;
    model.grid.dlon_deg = source_grid.dlon_deg;
    model.grid.dlat_deg = source_grid.dlat_deg;
    model.grid.nlon = nlon;
    model.grid.nlat = clipped_nlat;
    model.grid.global_lon = false;
    const std::size_t count = nlon * clipped_nlat;
    model.r2.resize(count);
    model.r1.resize(count);
    model.rho.resize(count);
    for (std::size_t i = 0; i < nlon; ++i) {
        const double lon = lon_first +
            static_cast<double>(i) * source_grid.dlon_deg;
        for (std::size_t j = 0; j < clipped_nlat; ++j) {
            const double lat = lat_first +
                static_cast<double>(j) * source_grid.dlat_deg;
            const State state = point_state(
                source_grid, source_r2, source_r1, source_rho, lon, lat);
            const std::size_t p = index(model.grid, i, j);
            model.r2[p] = state.r2;
            model.r1[p] = state.r1;
            model.rho[p] = state.rho;
        }
    }
    model.rectangle = {
        (observation_lon_deg - lon_half_width_deg) * degree,
        (observation_lon_deg + lon_half_width_deg) * degree,
        lat_south * degree, lat_north * degree};
    return model;
}

using HybridSurfaceKernel =
    double (*)(double, double, double, double, double, double);

inline Values hybrid_far_area(
    const Grid& g, const std::vector<double>& r2,
    const std::vector<double>& r1, const std::vector<double>& rho,
    const Point& point, const AngularRectangle& near_rectangle,
    double absolute_tolerance, double relative_tolerance,
    std::uint32_t canonical_component_mask) {
    static const std::array<HybridSurfaceKernel, 20> kernels{{
        surfacekernel_detail::KV, surfacekernel_detail::KVx,
        surfacekernel_detail::KVy, surfacekernel_detail::KVz,
        surfacekernel_detail::KVxx, surfacekernel_detail::KVxy,
        surfacekernel_detail::KVzx, surfacekernel_detail::KVyy,
        surfacekernel_detail::KVzy, surfacekernel_detail::KVzz,
        surfacekernel_detail::KVxxx, surfacekernel_detail::KVxxy,
        surfacekernel_detail::KVxxz, surfacekernel_detail::KVyyx,
        surfacekernel_detail::KVxyz, surfacekernel_detail::KVzzx,
        surfacekernel_detail::KVyyy, surfacekernel_detail::KVyyz,
        surfacekernel_detail::KVzzy, surfacekernel_detail::KVzzz}};
    const double observation_lat_deg = point.lat / degree;
    const double observation_lon_deg = point.lon / degree;
    const double west = near_rectangle.lon1 / degree;
    const double east = near_rectangle.lon2 / degree;
    const double south = near_rectangle.lat1 / degree;
    const double north = near_rectangle.lat2 / degree;
    const double domain_west = observation_lon_deg - 180.0;
    const double domain_east = observation_lon_deg + 180.0;
    const std::array<std::array<double, 4>, 4> rectangles{{
        {{-90.0, south, domain_west, domain_east}},
        {{north, 90.0, domain_west, domain_east}},
        {{south, north, domain_west, west}},
        {{south, north, east, domain_east}}}};
    Values result{};
    for (std::size_t k = 0; k < result.size(); ++k) {
        if ((canonical_component_mask & (1u << k)) == 0)
            continue;
        const HybridSurfaceKernel kernel = kernels[k];
        auto function = [&](const Matrix14& lat, const Matrix14& lon) {
            return surfacekernel_detail::EvaluateSurfaceKernel(
                r2, r1, rho, lat, lon, point.r,
                observation_lat_deg, observation_lon_deg,
                g.lon_first_deg, g.lon_last_deg,
                g.lat_first_deg, g.lat_last_deg,
                g.dlon_deg, g.dlat_deg, west, east, south, north,
                kernel);
        };
        double value = 0.0;
        for (const auto& rectangle : rectangles) {
            if (!(rectangle[1] > rectangle[0] &&
                  rectangle[3] > rectangle[2]))
                continue;
            Integral2Info info;
            value += integral2t(
                function, rectangle[0], rectangle[1],
                rectangle[2], rectangle[3], absolute_tolerance,
                relative_tolerance, 500, &info);
        }
        const std::size_t order =
            k == 0 ? 0 : (k <= 3 ? 1 : (k <= 9 ? 2 : 3));
        const std::array<double, 4> radial_scale{{
            point.r * point.r, point.r, 1.0, 1.0 / point.r}};
        result[k] = value * degree * degree * radial_scale[order];
    }
    return result;
}

inline Values point_estimate(
    const Grid& g, const std::vector<double>& r2,
    const std::vector<double>& r1, const std::vector<double>& rho,
    double observation_lat_deg, double observation_lon_deg,
    double observation_radius, double absolute_tolerance,
    double relative_tolerance, std::size_t component_context = 20,
    std::uint32_t canonical_component_mask = (1u << 20) - 1u) {
    const ScopedComponentContext requested_context(
        component_context, canonical_component_mask);
    const double normalized_lon_deg = detail::normalize_longitude(
        observation_lon_deg * degree) / degree;
    const Point point{observation_lat_deg * degree,
                      normalized_lon_deg * degree,
                      observation_radius};
    HybridNearModel near = make_hybrid_near_model(
        g, r2, r1, rho, normalized_lon_deg, observation_lat_deg);
    Values result{};
    if (near.grid.nlon == 1 && near.grid.nlat == 1) {
        const State state{near.r2[0], near.r1[0], near.rho[0]};
        result = constant_tesseroid_estimate(
            near.rectangle, state, point,
            absolute_tolerance, relative_tolerance);
    } else {
        result = point_estimate_tli(
            near.grid, near.r2, near.r1, near.rho,
            observation_lat_deg, normalized_lon_deg,
            observation_radius, absolute_tolerance, relative_tolerance,
            component_context, canonical_component_mask);
    }
    const Values far = hybrid_far_area(
        g, r2, r1, rho, point, near.rectangle,
        absolute_tolerance, relative_tolerance,
        canonical_component_mask);
    accumulate_values(far, result);
    return result;
}

inline void store(const Values& x, double* V, double* Vx, double* Vy,
                  double* Vz, double* Vxx, double* Vxy, double* Vyy,
                  double* Vzx, double* Vzy, double* Vzz, double* Vxxx,
                  double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx,
                  double* Vyyy, double* Vyyz, double* Vzzx, double* Vzzy,
                  double* Vzzz, std::uint32_t mask) {
    auto put = [&](std::uint32_t bit, double* output, std::size_t k) {
        if (output)
            *output = (mask & bit) ? x[k] : 0.0;
    };
    put(MSI_V, V, 0);
    put(MSI_Vx, Vx, 1);
    put(MSI_Vy, Vy, 2);
    put(MSI_Vz, Vz, 3);
    put(MSI_Vxx, Vxx, 4);
    put(MSI_Vxy, Vxy, 5);
    put(MSI_Vyy, Vyy, 7);
    put(MSI_Vzx, Vzx, 6);
    put(MSI_Vzy, Vzy, 8);
    put(MSI_Vzz, Vzz, 9);
    put(MSI_Vxxx, Vxxx, 10);
    put(MSI_Vxxy, Vxxy, 11);
    put(MSI_Vxxz, Vxxz, 12);
    put(MSI_Vxyz, Vxyz, 14);
    put(MSI_Vyyx, Vyyx, 13);
    put(MSI_Vyyy, Vyyy, 16);
    put(MSI_Vyyz, Vyyz, 17);
    put(MSI_Vzzx, Vzzx, 15);
    put(MSI_Vzzy, Vzzy, 18);
    put(MSI_Vzzz, Vzzz, 19);
}

inline std::size_t single_component_context(std::uint32_t mask) {
    if (mask == 0 || (mask & (mask - 1U)) != 0)
        return 20;
    static constexpr std::array<std::size_t, 20> canonical{{
        0, 1, 2, 3, 4, 5, 7, 6, 8, 9,
        10, 11, 12, 14, 13, 16, 17, 15, 18, 19}};
    std::size_t bit = 0;
    while ((mask & (1U << bit)) == 0)
        ++bit;
    return canonical[bit];
}

inline std::uint32_t canonical_component_mask(std::uint32_t public_mask) {
    static constexpr std::array<std::size_t, 20> canonical{{
        0, 1, 2, 3, 4, 5, 7, 6, 8, 9,
        10, 11, 12, 14, 13, 16, 17, 15, 18, 19}};
    std::uint32_t result = 0;
    for (std::size_t bit = 0; bit < canonical.size(); ++bit)
        if (public_mask & (1u << bit))
            result |= 1u << canonical[bit];
    return result;
}

}

inline void Tesseroid_GravityPointEstimate(
    double Lons, double Lone, double Lats, double Late, double dLon,
    double dLat, const std::vector<double>& R2,
    const std::vector<double>& R1, const std::vector<double>& rho, double fai,
    double lamda, double r, double AbsTol, double RelTol, double* V,
    double* Vx, double* Vy, double* Vz, double* Vxx, double* Vxy,
    double* Vyy, double* Vzx, double* Vzy, double* Vzz, double* Vxxx,
    double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx, double* Vyyy,
    double* Vyyz, double* Vzzx, double* Vzzy, double* Vzzz,
    std::uint32_t componentMask = MSI_All) {
    if (R2.size() != R1.size() || R2.size() != rho.size() || R2.empty())
        throw std::invalid_argument("TLI-MMSI source arrays must be non-empty and equally sized");
    const auto grid =
        tli_mmsi_detail::make_grid(Lons, Lone, Lats, Late, dLon, dLat, R2.size());
    const tli_mmsi_detail::Values values =
        tli_mmsi_detail::point_estimate(
            grid, R2, R1, rho, fai, lamda, r, AbsTol, RelTol,
            tli_mmsi_detail::single_component_context(componentMask),
            tli_mmsi_detail::canonical_component_mask(componentMask));
    tli_mmsi_detail::store(values, V, Vx, Vy, Vz, Vxx, Vxy, Vyy, Vzx,
                           Vzy, Vzz, Vxxx, Vxxy, Vxxz, Vxyz, Vyyx, Vyyy,
                           Vyyz, Vzzx, Vzzy, Vzzz, componentMask);
}

inline void Tesseroid_MagneticPointEstimate(
    double Lons, double Lone, double Lats, double Late, double dLon,
    double dLat, const std::vector<double>& R2,
    const std::vector<double>& R1, const std::vector<double>& Mx,
    const std::vector<double>& My, const std::vector<double>& Mz, double fai,
    double lamda, double r, double AbsTol, double RelTol, double* V,
    double* Vx, double* Vy, double* Vz, double* Vxx, double* Vxy,
    double* Vyy, double* Vzx, double* Vzy, double* Vzz) {
    if (Mx.size() != My.size() || Mx.size() != Mz.size() ||
        Mx.size() != R1.size() || Mx.size() != R2.size())
        throw std::invalid_argument("TLI-MMSI magnetization arrays must have matching sizes");

    const std::size_t nlat =
        static_cast<std::size_t>(std::llround((Late - Lats) / dLat)) + 1;
    const std::size_t nlon =
        static_cast<std::size_t>(std::llround((Lone - Lons) / dLon)) + 1;
    std::vector<double> mx_local(Mx.size()), my_local(My.size()),
        mz_local(Mz.size());
    const double co = std::cos(fai * tli_mmsi_detail::degree);
    const double so = std::sin(fai * tli_mmsi_detail::degree);
    for (std::size_t i = 0; i < nlon; ++i) {
        const double source_lon = Lons + static_cast<double>(i) * dLon;
        const double beta = (source_lon - lamda) * tli_mmsi_detail::degree;
        const double cb = std::cos(beta), sb = std::sin(beta);
        for (std::size_t j = 0; j < nlat; ++j) {
            const std::size_t p = i * nlat + j;
            const double source_lat =
                (Lats + static_cast<double>(j) * dLat) *
                tli_mmsi_detail::degree;
            const double cs = std::cos(source_lat), ss = std::sin(source_lat);
            mx_local[p] = (cb * so * ss + co * cs) * Mx[p] +
                          sb * so * My[p] -
                          (cb * so * cs - co * ss) * Mz[p];
            my_local[p] = -sb * ss * Mx[p] + cb * My[p] + sb * cs * Mz[p];
            mz_local[p] = -(cb * co * ss - so * cs) * Mx[p] -
                          sb * co * My[p] +
                          (cb * co * cs + so * ss) * Mz[p];
        }
    }

    std::array<double, 20> gx{}, gy{}, gz{};
    auto gravity = [&](const std::vector<double>& source,
                       std::array<double, 20>& out,
                       std::uint32_t component_mask) {
        Tesseroid_GravityPointEstimate(
            Lons, Lone, Lats, Late, dLon, dLat, R2, R1, source, fai, lamda,
            r, AbsTol, RelTol, &out[0], &out[1], &out[2], &out[3], &out[4],
            &out[5], &out[7], &out[6], &out[8], &out[9], &out[10], &out[11],
            &out[12], &out[14], &out[13], &out[16], &out[17], &out[15],
            &out[18], &out[19], component_mask);
    };
    constexpr std::uint32_t mx_mask =
        MSI_Vx | MSI_Vxx | MSI_Vxy | MSI_Vzx |
        MSI_Vxxx | MSI_Vxxy | MSI_Vyyx | MSI_Vxxz |
        MSI_Vxyz | MSI_Vzzx;
    constexpr std::uint32_t my_mask =
        MSI_Vy | MSI_Vxy | MSI_Vyy | MSI_Vzy |
        MSI_Vxxy | MSI_Vyyx | MSI_Vyyy | MSI_Vxyz |
        MSI_Vyyz | MSI_Vzzy;
    constexpr std::uint32_t mz_mask =
        MSI_Vz | MSI_Vzx | MSI_Vzy | MSI_Vzz |
        MSI_Vxxz | MSI_Vxyz | MSI_Vyyz | MSI_Vzzx |
        MSI_Vzzy | MSI_Vzzz;
    gravity(mx_local, gx, mx_mask);
    gravity(my_local, gy, my_mask);
    gravity(mz_local, gz, mz_mask);

    *V = gx[1] + gy[2] + gz[3];
    *Vx = gx[4] + gy[5] + gz[6];
    *Vy = gx[5] + gy[7] + gz[8];
    *Vz = gx[6] + gy[8] + gz[9];
    *Vxx = gx[10] + gy[11] + gz[12];
    *Vxy = gx[11] + gy[13] + gz[14];
    *Vyy = gx[13] + gy[16] + gz[17];
    *Vzx = gx[12] + gy[14] + gz[15];
    *Vzy = gx[14] + gy[17] + gz[18];
    *Vzz = gx[15] + gy[18] + gz[19];
}

inline void MassShellIntegration(
    const double Lons, const double Lone, const double Lats,
    const double Late, const double dLon, const double dLat,
    const double* R2In, const double* R1In, const double* rhoIn,
    const int NumPoint, const double* Latitude, const double* Longitude,
    const double* Radius, const double AbsTol, const double RelTol, double* V,
    double* Vx, double* Vy, double* Vz, double* Vxx, double* Vxy,
    double* Vyy, double* Vzx, double* Vzy, double* Vzz, double* Vxxx,
    double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx, double* Vyyy,
    double* Vyyz, double* Vzzx, double* Vzzy, double* Vzzz) {
    if (NumPoint <= 0)
        return;
    const std::size_t nlat =
        static_cast<std::size_t>(std::llround((Late - Lats) / dLat)) + 1;
    const std::size_t nlon =
        static_cast<std::size_t>(std::llround((Lone - Lons) / dLon)) + 1;
    const std::size_t count = nlon * nlat;
    const std::vector<double> r2(R2In, R2In + count);
    const std::vector<double> r1(R1In, R1In + count);
    const std::vector<double> rho(rhoIn, rhoIn + count);
    const auto grid = tli_mmsi_detail::make_grid(
        Lons, Lone, Lats, Late, dLon, dLat, count);
    for (int i = 0; i < NumPoint; ++i) {
        const auto values = tli_mmsi_detail::point_estimate(
            grid, r2, r1, rho, Latitude[i], Longitude[i], Radius[i], AbsTol,
            RelTol);
        tli_mmsi_detail::store(
            values, &V[i], &Vx[i], &Vy[i], &Vz[i], &Vxx[i], &Vxy[i],
            &Vyy[i], &Vzx[i], &Vzy[i], &Vzz[i], &Vxxx[i], &Vxxy[i],
            &Vxxz[i], &Vxyz[i], &Vyyx[i], &Vyyy[i], &Vyyz[i], &Vzzx[i],
            &Vzzy[i], &Vzzz[i], MSI_All);
    }
}

inline void MassShellIntegrationOpenMP(
    const double Lons, const double Lone, const double Lats,
    const double Late, const double dLon, const double dLat,
    const double* R2In, const double* R1In, const double* rhoIn,
    const int NumPoint, const double* Latitude, const double* Longitude,
    const double* Radius, const double AbsTol, const double RelTol, double* V,
    double* Vx, double* Vy, double* Vz, double* Vxx, double* Vxy,
    double* Vyy, double* Vzx, double* Vzy, double* Vzz, double* Vxxx,
    double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx, double* Vyyy,
    double* Vyyz, double* Vzzx, double* Vzzy, double* Vzzz) {
    if (NumPoint <= 0)
        return;
    const std::size_t nlat =
        static_cast<std::size_t>(std::llround((Late - Lats) / dLat)) + 1;
    const std::size_t nlon =
        static_cast<std::size_t>(std::llround((Lone - Lons) / dLon)) + 1;
    const std::size_t count = nlon * nlat;
    const std::vector<double> r2(R2In, R2In + count);
    const std::vector<double> r1(R1In, R1In + count);
    const std::vector<double> rho(rhoIn, rhoIn + count);
    const auto grid = tli_mmsi_detail::make_grid(
        Lons, Lone, Lats, Late, dLon, dLat, count);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(NumPoint > 1)
#endif
    for (int i = 0; i < NumPoint; ++i) {
        const auto values = tli_mmsi_detail::point_estimate(
            grid, r2, r1, rho, Latitude[i], Longitude[i], Radius[i], AbsTol,
            RelTol);
        tli_mmsi_detail::store(
            values, &V[i], &Vx[i], &Vy[i], &Vz[i], &Vxx[i], &Vxy[i],
            &Vyy[i], &Vzx[i], &Vzy[i], &Vzz[i], &Vxxx[i], &Vxxy[i],
            &Vxxz[i], &Vxyz[i], &Vyyx[i], &Vyyy[i], &Vyyz[i], &Vzzx[i],
            &Vzzy[i], &Vzzz[i], MSI_All);
    }
}

inline void MagnetizationShellIntegration(
    const double Lons, const double Lone, const double Lats,
    const double Late, const double dLon, const double dLat,
    const double* R2In, const double* R1In, const double* MxIn,
    const double* MyIn, const double* MzIn, const int NumPoint,
    const double* Latitude, const double* Longitude, const double* Radius,
    const double AbsTol, const double RelTol, double* V, double* Vx,
    double* Vy, double* Vz, double* Vxx, double* Vxy, double* Vyy,
    double* Vzx, double* Vzy, double* Vzz) {
    if (NumPoint <= 0)
        return;
    const std::size_t nlat =
        static_cast<std::size_t>(std::llround((Late - Lats) / dLat)) + 1;
    const std::size_t nlon =
        static_cast<std::size_t>(std::llround((Lone - Lons) / dLon)) + 1;
    const std::size_t count = nlon * nlat;
    const std::vector<double> r2(R2In, R2In + count);
    const std::vector<double> r1(R1In, R1In + count);
    const std::vector<double> mx(MxIn, MxIn + count);
    const std::vector<double> my(MyIn, MyIn + count);
    const std::vector<double> mz(MzIn, MzIn + count);
    for (int i = 0; i < NumPoint; ++i)
        Tesseroid_MagneticPointEstimate(
            Lons, Lone, Lats, Late, dLon, dLat, r2, r1, mx, my, mz,
            Latitude[i], Longitude[i], Radius[i], AbsTol, RelTol, &V[i],
            &Vx[i], &Vy[i], &Vz[i], &Vxx[i], &Vxy[i], &Vyy[i], &Vzx[i],
            &Vzy[i], &Vzz[i]);
}

inline void MagnetizationShellIntegrationOpenMP(
    const double Lons, const double Lone, const double Lats,
    const double Late, const double dLon, const double dLat,
    const double* R2In, const double* R1In, const double* MxIn,
    const double* MyIn, const double* MzIn, const int NumPoint,
    const double* Latitude, const double* Longitude, const double* Radius,
    const double AbsTol, const double RelTol, double* V, double* Vx,
    double* Vy, double* Vz, double* Vxx, double* Vxy, double* Vyy,
    double* Vzx, double* Vzy, double* Vzz) {
    if (NumPoint <= 0)
        return;
    const std::size_t nlat =
        static_cast<std::size_t>(std::llround((Late - Lats) / dLat)) + 1;
    const std::size_t nlon =
        static_cast<std::size_t>(std::llround((Lone - Lons) / dLon)) + 1;
    const std::size_t count = nlon * nlat;
    const std::vector<double> r2(R2In, R2In + count);
    const std::vector<double> r1(R1In, R1In + count);
    const std::vector<double> mx(MxIn, MxIn + count);
    const std::vector<double> my(MyIn, MyIn + count);
    const std::vector<double> mz(MzIn, MzIn + count);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < NumPoint; ++i)
        Tesseroid_MagneticPointEstimate(
            Lons, Lone, Lats, Late, dLon, dLat, r2, r1, mx, my, mz,
            Latitude[i], Longitude[i], Radius[i], AbsTol, RelTol, &V[i],
            &Vx[i], &Vy[i], &Vz[i], &Vxx[i], &Vxy[i], &Vyy[i], &Vzx[i],
            &Vzy[i], &Vzz[i]);
}

#ifdef TLI_ENABLE_MPI
inline void MassShellIntegrationMPI(
    const double Lons, const double Lone, const double Lats,
    const double Late, const double dLon, const double dLat,
    const double* R2In, const double* R1In, const double* rhoIn,
    const int NumPoint, const double* Latitude, const double* Longitude,
    const double* Radius, const double AbsTol, const double RelTol, double* V,
    double* Vx, double* Vy, double* Vz, double* Vxx, double* Vxy,
    double* Vyy, double* Vzx, double* Vzy, double* Vzz, double* Vxxx,
    double* Vxxy, double* Vxxz, double* Vxyz, double* Vyyx, double* Vyyy,
    double* Vyyz, double* Vzzx, double* Vzzy, double* Vzzz) {
    int rank = 0, processes = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    const std::size_t nlat =
        static_cast<std::size_t>(std::llround((Late - Lats) / dLat)) + 1;
    const std::size_t nlon =
        static_cast<std::size_t>(std::llround((Lone - Lons) / dLon)) + 1;
    const std::size_t count = nlon * nlat;
    std::vector<double> r2(count), r1(count), rho(count);
    std::vector<double> latitude(NumPoint), longitude(NumPoint), radius(NumPoint);
    if (rank == 0) {
        std::copy(R2In, R2In + count, r2.begin());
        std::copy(R1In, R1In + count, r1.begin());
        std::copy(rhoIn, rhoIn + count, rho.begin());
        std::copy(Latitude, Latitude + NumPoint, latitude.begin());
        std::copy(Longitude, Longitude + NumPoint, longitude.begin());
        std::copy(Radius, Radius + NumPoint, radius.begin());
    }
    MPI_Bcast(r2.data(), static_cast<int>(count), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(r1.data(), static_cast<int>(count), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rho.data(), static_cast<int>(count), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(latitude.data(), NumPoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(longitude.data(), NumPoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(radius.data(), NumPoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::vector<int> counts(processes), offsets(processes);
    for (int p = 0; p < processes; ++p) {
        counts[p] = NumPoint / processes + (p < NumPoint % processes ? 1 : 0);
        offsets[p] = p == 0 ? 0 : offsets[p - 1] + counts[p - 1];
    }
    const int local_count = counts[rank], offset = offsets[rank];
    std::array<std::vector<double>, 20> local;
    for (auto& values : local)
        values.resize(local_count);
    for (int i = 0; i < local_count; ++i) {
        const int p = offset + i;
        Tesseroid_GravityPointEstimate(
            Lons, Lone, Lats, Late, dLon, dLat, r2, r1, rho, latitude[p],
            longitude[p], radius[p], AbsTol, RelTol, &local[0][i],
            &local[1][i], &local[2][i], &local[3][i], &local[4][i],
            &local[5][i], &local[7][i], &local[6][i], &local[8][i],
            &local[9][i], &local[10][i], &local[11][i], &local[12][i],
            &local[14][i], &local[13][i], &local[16][i], &local[17][i],
            &local[15][i], &local[18][i], &local[19][i]);
    }
    const std::array<double*, 20> output{
        V, Vx, Vy, Vz, Vxx, Vxy, Vzx, Vyy, Vzy, Vzz,
        Vxxx, Vxxy, Vxxz, Vyyx, Vxyz, Vzzx, Vyyy, Vyyz, Vzzy, Vzzz};
    for (std::size_t k = 0; k < output.size(); ++k)
        MPI_Gatherv(local[k].data(), local_count, MPI_DOUBLE, output[k],
                    counts.data(), offsets.data(), MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);
}

inline void MagnetizationShellIntegrationMPI(
    const double Lons, const double Lone, const double Lats,
    const double Late, const double dLon, const double dLat,
    const double* R2In, const double* R1In, const double* MxIn,
    const double* MyIn, const double* MzIn, const int NumPoint,
    const double* Latitude, const double* Longitude, const double* Radius,
    const double AbsTol, const double RelTol, double* V, double* Vx,
    double* Vy, double* Vz, double* Vxx, double* Vxy, double* Vyy,
    double* Vzx, double* Vzy, double* Vzz) {
    int rank = 0, processes = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    const std::size_t nlat =
        static_cast<std::size_t>(std::llround((Late - Lats) / dLat)) + 1;
    const std::size_t nlon =
        static_cast<std::size_t>(std::llround((Lone - Lons) / dLon)) + 1;
    const std::size_t count = nlon * nlat;
    std::vector<double> r2(count), r1(count), mx(count), my(count), mz(count);
    std::vector<double> latitude(NumPoint), longitude(NumPoint), radius(NumPoint);
    if (rank == 0) {
        std::copy(R2In, R2In + count, r2.begin());
        std::copy(R1In, R1In + count, r1.begin());
        std::copy(MxIn, MxIn + count, mx.begin());
        std::copy(MyIn, MyIn + count, my.begin());
        std::copy(MzIn, MzIn + count, mz.begin());
        std::copy(Latitude, Latitude + NumPoint, latitude.begin());
        std::copy(Longitude, Longitude + NumPoint, longitude.begin());
        std::copy(Radius, Radius + NumPoint, radius.begin());
    }
    for (double* data : {r2.data(), r1.data(), mx.data(), my.data(), mz.data()})
        MPI_Bcast(data, static_cast<int>(count), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(latitude.data(), NumPoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(longitude.data(), NumPoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(radius.data(), NumPoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::vector<int> counts(processes), offsets(processes);
    for (int p = 0; p < processes; ++p) {
        counts[p] = NumPoint / processes + (p < NumPoint % processes ? 1 : 0);
        offsets[p] = p == 0 ? 0 : offsets[p - 1] + counts[p - 1];
    }
    const int local_count = counts[rank], offset = offsets[rank];
    std::array<std::vector<double>, 10> local;
    for (auto& values : local)
        values.resize(local_count);
    for (int i = 0; i < local_count; ++i) {
        const int p = offset + i;
        Tesseroid_MagneticPointEstimate(
            Lons, Lone, Lats, Late, dLon, dLat, r2, r1, mx, my, mz,
            latitude[p], longitude[p], radius[p], AbsTol, RelTol,
            &local[0][i], &local[1][i], &local[2][i], &local[3][i],
            &local[4][i], &local[5][i], &local[6][i], &local[7][i],
            &local[8][i], &local[9][i]);
    }
    const std::array<double*, 10> output{
        V, Vx, Vy, Vz, Vxx, Vxy, Vyy, Vzx, Vzy, Vzz};
    for (std::size_t k = 0; k < output.size(); ++k)
        MPI_Gatherv(local[k].data(), local_count, MPI_DOUBLE, output[k],
                    counts.data(), offsets.data(), MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);
}
#endif

}
