#include "../TLI_MMSI/_tesseroid_MSI.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
    using namespace Tesseroid;
    constexpr double degree = PI / 180.0;

    EstimateOptions options;
    options.integration.atol = 1.0e-14;
    options.integration.rtol = 1.0e-9;
    options.gravitational_constant = 1.0;
    const Prism prism{34.95 * degree, 35.05 * degree,
                      119.95 * degree, 120.05 * degree,
                      6341.2, 6371.2, 2670.0};
    const Point point{35.051 * degree, 120.0 * degree, 6372.2};
    const auto tli = estimate_all(prism, point, options);
    for (const auto& component : tli)
        if (!std::isfinite(component.value) || !std::isfinite(component.error))
            return 1;

    const auto grid = tli_mmsi_detail::make_grid(
        119.95, 120.05, 34.95, 35.05, 0.1, 0.1, 4);
    const std::vector<double> r2(4, 6371.2);
    const std::vector<double> r1(4, 6341.2);
    const std::vector<double> density(4, 2670.0);
    const auto mmsi = tli_mmsi_detail::point_estimate(
        grid, r2, r1, density, 35.0, 120.0, 6372.2,
        1.0e-14, 1.0e-9);
    for (double value : mmsi)
        if (!std::isfinite(value))
            return 2;

    const double trace2 = mmsi[4] + mmsi[7] + mmsi[9];
    const double scale2 = std::max({1.0, std::abs(mmsi[4]),
                                    std::abs(mmsi[7]), std::abs(mmsi[9])});
    if (std::abs(trace2) / scale2 > 1.0e-8)
        return 3;

    std::cout << "TLI and TLI-MMSI smoke test passed\n";
    return 0;
}
