// Small-scale timing probe for the TLI line-integration kernel.
// Measures single-component wall-clock times at a few representative
// observation points of the RunTimeTest1 benchmark grid, so the full
// 160,801-point cluster run can be extrapolated without launching MPI.
//
// Points probed (all at RE + 1 km):
//   1. far outside the prism  (27.0N, 105.0E)
//   2. near angular boundary  (30.2N, 110.2E, ~0.15 deg from corner)
//   3. inside the projection  (40.0N, 120.0E)
//   4. near the prism centre with a displaced longitude (40.0N, 129.5E)
//   5. near the other corner  (49.8N, 129.8E)
//
// Build (from the RunTimeTest1 directory):
//   icpx -O2 -std=c++17 -I. Timing_Probe.cpp -o Timing_Probe
// Run:
//   Timing_Probe.exe

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "_TLI_estimate.h"
#include "_TLI_integralkernel.h"
#include "_integral.h"

namespace
{

constexpr double EARTH_RADIUS_M = 6371200.0;
constexpr double OBSERVATION_HEIGHT_M = 1000.0;
constexpr double PRISM_LAT_MIN_DEG = 30.0;
constexpr double PRISM_LAT_MAX_DEG = 50.0;
constexpr double PRISM_LON_MIN_DEG = 110.0;
constexpr double PRISM_LON_MAX_DEG = 130.0;
constexpr double PRISM_R_MIN_M = EARTH_RADIUS_M - 30000.0;
constexpr double PRISM_R_MAX_M = EARTH_RADIUS_M;
constexpr double PRISM_DENSITY = 2670.0;

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double DEG_TO_RAD = PI / 180.0;

const std::vector<std::pair<double, double>> PROBE_POINTS = {
    {27.0, 105.0},   // far outside
    {30.2, 110.2},   // near boundary
    {40.0, 120.0},   // inside
    {40.0, 129.5},   // inside, near east boundary
    {49.8, 129.8}    // near corner
};

const std::vector<const char*> COMPONENT_NAMES = {
    "V", "Vx", "Vy", "Vz",
    "Vxx", "Vxy", "Vyy", "Vzx", "Vzy", "Vzz",
    "Vxxx", "Vxxy", "Vxxz", "Vxyz", "Vyyx",
    "Vyyy", "Vyyz", "Vzzx", "Vzzy", "Vzzz"
};

using Clock = std::chrono::steady_clock;

double evaluateOneComponent(
    const char* componentName, double latitude, double longitude,
    double radius, bool& converged)
{
    using namespace Tesseroid;

    const Prism prism{
        PRISM_LAT_MIN_DEG * DEG_TO_RAD,
        PRISM_LAT_MAX_DEG * DEG_TO_RAD,
        PRISM_LON_MIN_DEG * DEG_TO_RAD,
        PRISM_LON_MAX_DEG * DEG_TO_RAD,
        PRISM_R_MIN_M,
        PRISM_R_MAX_M,
        PRISM_DENSITY};

    const Point point{
        latitude * DEG_TO_RAD,
        detail::normalize_longitude(longitude * DEG_TO_RAD),
        radius};

    EstimateOptions options;
    options.integration.atol = 1.0e-16;
    options.integration.rtol = 1.0e-16;
    options.gravitational_constant = 1.0;

    EstimateInfo info;
    const double value = estimate(component(componentName), prism, point,
                                  options, &info);
    converged = info.converged;
    return value;
}

} // namespace

int main()
{
    const double radius = EARTH_RADIUS_M + OBSERVATION_HEIGHT_M;

    std::cout << std::setprecision(6)
              << "point,lat,lon,component,mean_ms,value,converged\n";

    double worstComponentMs = 0.0;
    double totalPointMs = 0.0;
    std::vector<double> componentMaxMs(20, 0.0);

    for (std::size_t p = 0; p < PROBE_POINTS.size(); ++p) {
        const double lat = PROBE_POINTS[p].first;
        const double lon = PROBE_POINTS[p].second;
        double pointMs = 0.0;

        for (std::size_t c = 0; c < COMPONENT_NAMES.size(); ++c) {
            // One warm-up evaluation keeps instruction and data caches hot.
            bool warmConverged = false;
            evaluateOneComponent(COMPONENT_NAMES[c], lat, lon, radius,
                                 warmConverged);

            // Two timed evaluations for a stable mean.
            double bestMs = 1.0e30;
            double value = 0.0;
            bool converged = false;
            for (int repeat = 0; repeat < 2; ++repeat) {
                const auto start = Clock::now();
                value = evaluateOneComponent(COMPONENT_NAMES[c], lat, lon,
                                             radius, converged);
                const double ms =
                    std::chrono::duration<double, std::milli>(
                        Clock::now() - start).count();
                if (ms < bestMs) bestMs = ms;
            }

            pointMs += bestMs;
            componentMaxMs[c] = std::max(componentMaxMs[c], bestMs);
            worstComponentMs = std::max(worstComponentMs, bestMs);

            std::cout << p + 1 << ',' << lat << ',' << lon << ','
                      << COMPONENT_NAMES[c] << ',' << bestMs << ','
                      << value << ',' << (converged ? 1 : 0) << '\n';
        }

        std::cout << p + 1 << ',' << lat << ',' << lon
                  << ",POINT_TOTAL," << pointMs << ",,\n";
        totalPointMs += pointMs;
    }

    const double meanPointMs = totalPointMs / PROBE_POINTS.size();
    std::cout << "\n=== Summary (best of 2 timed runs, warmup included) ===\n"
              << "Mean 20-component point time : " << meanPointMs << " ms\n"
              << "Worst single component time  : " << worstComponentMs << " ms\n";

    std::cout << "Per-component maxima (ms):\n";
    for (std::size_t c = 0; c < COMPONENT_NAMES.size(); ++c)
        std::cout << "  " << COMPONENT_NAMES[c] << ": "
                  << componentMaxMs[c] << '\n';

    // Extrapolation for the full benchmark:
    //   160,801 points x (2 warmup + 10 timed) evaluations per component
    //   x 20 components per point, on one core.
    constexpr double GRID_POINTS = 160801.0;
    constexpr double EVALUATIONS_PER_COMPONENT = 12.0;
    const double singleCoreSeconds =
        meanPointMs * 1.0e-3 * EVALUATIONS_PER_COMPONENT * GRID_POINTS;
    std::cout << "\n=== Extrapolation (single core, 160,801 points, 12x per component) ===\n"
              << "Single-core CPU time : " << singleCoreSeconds / 3600.0
              << " h = " << singleCoreSeconds / 86400.0 << " days\n";
    for (double processes : {1800.0, 900.0, 180.0, 45.0}) {
        std::cout << "  ~" << static_cast<int>(processes)
                  << " MPI processes wall-clock (perfect load balance): "
                  << singleCoreSeconds / processes / 3600.0 << " h\n";
    }

    return 0;
}
