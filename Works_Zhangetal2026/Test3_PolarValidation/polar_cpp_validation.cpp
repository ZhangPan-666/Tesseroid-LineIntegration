// Test 3 numerical side: ordinary Tesseroid C++ implementation.
#include "_TLI_estimate.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

int main() {
    constexpr int n = 2001;
    constexpr double re = 6371.2;
    const std::vector<std::string> names{"Normal", "Small", "Thin", "ZonalBand"};
    const std::vector<std::array<double, 6>> p{
        {69.99, 89.99, -10.0, 10.0, 6341.2, 6371.2},
        {89.89, 89.99, -0.05, 0.05, 6370.2, 6371.2},
        {69.99, 89.99, -10.0, 10.0, 6371.19, 6371.2},
        {69.99, 89.99, -180.0, 180.0, 6341.2, 6371.2}};
    std::vector<double> fai1(4), fai2(4), lon1(4), lon2(4), r1(4), r2(4), rho(4, 1.0);
    for (int k = 0; k < 4; ++k) {
        fai1[k] = p[k][0]; fai2[k] = p[k][1]; lon1[k] = p[k][2];
        lon2[k] = p[k][3]; r1[k] = p[k][4]; r2[k] = p[k][5];
    }
    std::ofstream out("results/polar_cpp_results.csv");
    if (!out) return 2;
    out << "case,height,V,Vz,Vxx,Vyy,Vzz,Vxxz,Vyyz,Vzzz\n";
    out << std::setprecision(17) << std::scientific;
    std::vector<double> lat(n, 90.0), lon(n, 0.0), radius(n);
    std::vector<double> V(n), Vx(n), Vy(n), Vz(n), Vxx(n), Vxy(n), Vyy(n), Vzx(n), Vzy(n), Vzz(n);
    std::vector<double> Vxxx(n), Vxxy(n), Vxxz(n), Vxyz(n), Vyyx(n), Vyyy(n), Vyyz(n), Vzzx(n), Vzzy(n), Vzzz(n);
    for (int i = 0; i < n; ++i) radius[i] = re + static_cast<double>(i - 1000);
    for (int k = 0; k < 4; ++k) {
        Tesseroid::Tesseroid_GravityEstimate(
            1, &fai1[k], &fai2[k], &lon1[k], &lon2[k], &r1[k], &r2[k], &rho[k], n,
            lat.data(), lon.data(), radius.data(), 1e-14, 1e-14,
            V.data(), Vx.data(), Vy.data(), Vz.data(), Vxx.data(), Vxy.data(), Vyy.data(),
            Vzx.data(), Vzy.data(), Vzz.data(), Vxxx.data(), Vxxy.data(), Vxxz.data(),
            Vxyz.data(), Vyyx.data(), Vyyy.data(), Vyyz.data(), Vzzx.data(), Vzzy.data(), Vzzz.data());
        for (int i = 0; i < n; ++i)
            out << names[k] << ',' << (i - 1000) << ',' << V[i] << ',' << Vz[i] << ','
                << Vxx[i] << ',' << Vyy[i] << ',' << Vzz[i] << ',' << Vxxz[i] << ','
                << Vyyz[i] << ',' << Vzzz[i] << '\n';
    }
    return 0;
}
