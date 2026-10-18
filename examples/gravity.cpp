#include "_TLI_estimate.h"
#include <iomanip>
#include <iostream>

int main() {
    using namespace Tesseroid;
    constexpr double deg = PI / 180.0;
    // SI: metres, kg/m^3, and the default G in m^3 kg^-1 s^-2.
    const Prism source{34.8*deg, 35.2*deg, 119.7*deg, 120.3*deg,
                       6341200.0, 6371200.0, 2670.0};
    const Point observation{35.11*deg, 120.08*deg, 6771200.0};
    EstimateOptions options;
    options.integration.atol = 1e-18;
    options.integration.rtol = 1e-9;
    const auto result = estimate_all(source, observation, options);
    std::cout << "component,value,estimated_absolute_error,converged\n"
              << std::scientific << std::setprecision(15);
    bool converged = true;
    for (const auto& item : result) {
        std::cout << component_name(item.component) << ',' << item.value << ','
                  << item.error << ',' << item.converged << '\n';
        converged = converged && item.converged;
    }
    return converged ? 0 : 1;
}
