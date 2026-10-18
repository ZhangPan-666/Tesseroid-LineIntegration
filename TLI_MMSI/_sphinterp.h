#pragma once
#include "_integral2.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace sphinterp {

template <typename T> inline T clamp_value(const T& x, const T& lo, const T& hi) {
    return std::max(lo, std::min(x, hi));
}

inline double NormalizeLon180(double lon) {
    while (lon < -180.0)
        lon += 360.0;
    while (lon >= 180.0)
        lon -= 360.0;
    return lon;
}

inline bool InWrappedLonInterval(double lon, double Lon0s, double Lon0e) {
    constexpr double epsLon = 1e-12;

    if ((std::abs(Lon0s + 180.0) <= epsLon && std::abs(Lon0e - 180.0) <= epsLon) ||
        (std::abs(Lon0s - Lon0e) >= 360.0 - epsLon)) {
        return true;
    }

    lon = NormalizeLon180(lon);
    Lon0s = NormalizeLon180(Lon0s);
    Lon0e = NormalizeLon180(Lon0e);

    if (Lon0s <= Lon0e) {
        return (lon >= Lon0s && lon <= Lon0e);
    } else {
        return (lon >= Lon0s || lon <= Lon0e);
    }
}

template <std::size_t R = 14, std::size_t C = 14>
inline void SphBilinearInterp(const double Lons, const double Lone, const double Lats, const double Late,
                              const double dLon, const double dLat, const std::vector<double>& F1,
                              const std::vector<double>& F2, const std::vector<double>& F3,
                              const Tesseroid::FixedMatrix<R, C>& LonQ,
                              const Tesseroid::FixedMatrix<R, C>& LatQ,
                              const double Lon0s, const double Lon0e, const double Lat0s, const double Lat0e,
                              Tesseroid::FixedMatrix<R, C>& O1,
                              Tesseroid::FixedMatrix<R, C>& O2,
                              Tesseroid::FixedMatrix<R, C>& O3) {
    const long long NumLat = static_cast<long long>(std::llround((Late - Lats) / dLat)) + 1;
    const long long NumLon = static_cast<long long>(std::llround((Lone - Lons) / dLon)) + 1;

    O1.fill(0.0);
    O2.fill(0.0);
    O3.fill(0.0);

    const double pi = 4.0 * std::atan(1.0);
    const double deg2rad = pi / 180.0;

    for (std::size_t i = 0; i < R; ++i) {
        for (std::size_t j = 0; j < C; ++j) {

            double lon0 = NormalizeLon180(LonQ(i, j));
            const double lat0 = LatQ(i, j);

            if (InWrappedLonInterval(lon0, Lon0s, Lon0e) && lat0 >= Lat0s && lat0 <= Lat0e) {
                O1(i, j) = 0.0;
                O2(i, j) = 0.0;
                O3(i, j) = 0.0;
                continue;
            }

            if (lon0 < Lons || lon0 > Lone || lat0 < Lats || lat0 > Late) {

                const long long k0 =
                    clamp_value(static_cast<long long>(std::llround((lon0 - Lons) / dLon)) + 1, 1LL, NumLon);

                const long long l0 =
                    clamp_value(static_cast<long long>(std::llround((lat0 - Lats) / dLat)) + 1, 1LL, NumLat);

                const long long idx = (k0 - 1) * NumLat + l0 - 1;

                O1(i, j) = F1[idx];
                O2(i, j) = F2[idx];
                O3(i, j) = F3[idx];
                continue;
            }

            const double epsIndex = 1e-12;

            long long l = static_cast<long long>(std::floor((lat0 - Lats) / dLat + epsIndex)) + 1;
            long long k = static_cast<long long>(std::floor((lon0 - Lons) / dLon + epsIndex)) + 1;

            l = clamp_value(l, 1LL, NumLat - 1);
            k = clamp_value(k, 1LL, NumLon - 1);

            const double lat1 = Lats + (l - 1) * dLat;
            const double lat2 = lat1 + dLat;

            const double x1 = std::cos((90.0 - lat1) * deg2rad);
            const double x2 = std::cos((90.0 - lat2) * deg2rad);

            const double lon1 = Lons + (k - 1) * dLon;
            const double lon2 = lon1 + dLon;

            const double FF1 = lon1 * deg2rad;
            const double FF2 = lon2 * deg2rad;

            const double x0 = std::cos((90.0 - lat0) * deg2rad);
            const double phi0 = lon0 * deg2rad;

            const long long idx00 = (k - 1) * NumLat + l - 1;
            const long long idx10 = k * NumLat + l - 1;
            const long long idx01 = (k - 1) * NumLat + l;
            const long long idx11 = k * NumLat + l;

            const double dx = x2 - x1;
            const double dphi = FF2 - FF1;
            const double dd = dx * dphi;

            const double epsDD = 1e-14;

            auto interpOneField = [&](const std::vector<double>& F) -> double {
                const double f00 = F[idx00];
                const double f10 = F[idx10];
                const double f01 = F[idx01];
                const double f11 = F[idx11];

                if (std::abs(dx) < epsDD && std::abs(dphi) < epsDD) {
                    return f00;
                }

                if (std::abs(dx) < epsDD) {
                    const double t = (phi0 - FF1) / dphi;
                    return f00 + t * (f10 - f00);
                }

                if (std::abs(dphi) < epsDD) {
                    const double t = (x0 - x1) / dx;
                    return f00 + t * (f01 - f00);
                }

                if (std::abs(dd) < epsDD) {
                    const double tphi = (phi0 - FF1) / dphi;
                    const double tx = (x0 - x1) / dx;
                    const bool right = (tphi >= 0.5);
                    const bool top = (tx >= 0.5);
                    if (!right && !top)
                        return f00;
                    if (right && !top)
                        return f10;
                    if (!right && top)
                        return f01;
                    return f11;
                }

                const double tphi = (phi0 - FF1) / dphi;
                const double tx = (x0 - x1) / dx;
                const double bottom = f00 + tphi * (f10 - f00);
                const double top = f01 + tphi * (f11 - f01);
                return bottom + tx * (top - bottom);
            };

            O1(i, j) = interpOneField(F1);
            O2(i, j) = interpOneField(F2);
            O3(i, j) = interpOneField(F3);
        }
    }
}

}
