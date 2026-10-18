#pragma once
#include "_sphinterp.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace Phithreshold {

const double PI = 4.0 * std::atan(1.0);
const double Deg2Rad = PI / 180.0;
}

namespace Tesseroid {

namespace surfacekernel_detail {
inline double ClampCos(const double x) { return std::max(-1.0, std::min(1.0, x)); }

inline double Binom(const int n, const int k) {
    if (k < 0 || k > n)
        return 0.0;
    if (k == 0 || k == n)
        return 1.0;

    const int kk = std::min(k, n - k);
    double out = 1.0;
    for (int i = 1; i <= kk; ++i) {
        out *= static_cast<double>(n - kk + i) / static_cast<double>(i);
    }
    return out;
}

inline void LocalGeometry(const double FaiDeg, const double LamdaDeg, const double faiDeg, const double lamdaDeg,
                          double& Tx, double& Ty, double& SinPsi, double& CosPsi) {
    const double fai = faiDeg * Phithreshold::Deg2Rad;
    const double Fai = FaiDeg * Phithreshold::Deg2Rad;
    const double dLam = (LamdaDeg - lamdaDeg) * Phithreshold::Deg2Rad;

    const double cos_fai = std::cos(fai);
    const double sin_fai = std::sin(fai);
    const double cos_Fai = std::cos(Fai);
    const double sin_Fai = std::sin(Fai);
    const double cos_dLam = std::cos(dLam);
    const double sin_dLam = std::sin(dLam);

    Tx = cos_fai * sin_Fai - sin_fai * cos_Fai * cos_dLam;
    Ty = cos_Fai * sin_dLam;

    SinPsi = std::hypot(Tx, Ty);
    CosPsi = sin_fai * sin_Fai + cos_fai * cos_Fai * cos_dLam;
    CosPsi = ClampCos(CosPsi);
}

inline double LogRatio(const double a, const double b) {
    const double r = (a - b) / b;
    if (std::abs(r) < 0.25) {
        return std::log1p(r);
    }
    return std::log(a / b);
}

inline double StablePowerIntegral(const int alpha, const double q1, const double q2) {
    if (q1 <= 0.0 || q2 <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (q1 == q2) {
        return 0.0;
    }

    const double logRatio = LogRatio(q2, q1);

    if (alpha == -1) {
        return logRatio;
    }

    const int betaInt = alpha + 1;
    const double beta = static_cast<double>(betaInt);
    return std::pow(q1, beta) * std::expm1(beta * logRatio) / beta;
}

inline double BaseMomentSmallB(const double u1, const double u2, const double b, const int n, const int k) {
    if (u1 * u2 <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double sigma = (u1 > 0.0) ? 1.0 : -1.0;
    const double q1 = std::min(std::abs(u1), std::abs(u2));
    const double q2 = std::max(std::abs(u1), std::abs(u2));

    const double P0 = StablePowerIntegral(k - n, q1, q2);
    const double P2 = StablePowerIntegral(k - n - 2, q1, q2);
    const double P4 = StablePowerIntegral(k - n - 4, q1, q2);
    const double P6 = StablePowerIntegral(k - n - 6, q1, q2);
    const double P8 = StablePowerIntegral(k - n - 8, q1, q2);

    const double signFactor = (k % 2 == 0) ? 1.0 : sigma;
    return signFactor * (P0 - 0.5 * n * b * b * P2 + n * (n + 2.0) / 8.0 * std::pow(b, 4) * P4 -
                         n * (n + 2.0) * (n + 4.0) / 48.0 * std::pow(b, 6) * P6 +
                         n * (n + 2.0) * (n + 4.0) * (n + 6.0) / 384.0 * std::pow(b, 8) * P8);
}

inline double BaseIntegral(const double t, const double b, const int n, const int k) {
    const double L = std::sqrt(t * t + b * b);
    const double A = std::asinh(t / b);

    switch (n) {
    case 1:
        switch (k) {
        case 0:
            return A;
        case 1:
            return L;
        case 2:
            return 0.5 * (t * L - b * b * A);
        default:
            throw std::runtime_error("Unsupported k for n=1.");
        }

    case 3: {
        const double I30 = t / (b * b * L);
        const double I31 = -1.0 / L;
        const double I32 = A - b * b * I30;
        const double I33 = L + b * b / L;

        switch (k) {
        case 0:
            return I30;
        case 1:
            return I31;
        case 2:
            return I32;
        case 3:
            return I33;
        default:
            throw std::runtime_error("Unsupported k for n=3.");
        }
    }

    case 5: {
        const double L3 = L * L * L;
        const double I30 = t / (b * b * L);
        const double I31 = -1.0 / L;
        const double I32 = A - b * b * I30;

        const double I50 = t / (3.0 * b * b * L3) + 2.0 / (3.0 * b * b) * I30;
        const double I51 = -1.0 / (3.0 * L3);
        const double I52 = I30 - b * b * I50;
        const double I53 = I31 - b * b * I51;
        const double I54 = I32 - b * b * I52;

        switch (k) {
        case 0:
            return I50;
        case 1:
            return I51;
        case 2:
            return I52;
        case 3:
            return I53;
        case 4:
            return I54;
        default:
            throw std::runtime_error("Unsupported k for n=5.");
        }
    }

    case 7: {
        const double L3 = L * L * L;
        const double L5 = L3 * L * L;

        const double I30 = t / (b * b * L);
        const double I31 = -1.0 / L;
        const double I32 = A - b * b * I30;

        const double I50 = t / (3.0 * b * b * L3) + 2.0 / (3.0 * b * b) * I30;
        const double I51 = -1.0 / (3.0 * L3);
        const double I52 = I30 - b * b * I50;
        const double I53 = I31 - b * b * I51;
        const double I54 = I32 - b * b * I52;

        const double I70 = t / (5.0 * b * b * L5) + 4.0 / (5.0 * b * b) * I50;
        const double I71 = -1.0 / (5.0 * L5);
        const double I72 = I50 - b * b * I70;
        const double I73 = I51 - b * b * I71;
        const double I74 = I52 - b * b * I72;
        const double I75 = I53 - b * b * I73;

        switch (k) {
        case 0:
            return I70;
        case 1:
            return I71;
        case 2:
            return I72;
        case 3:
            return I73;
        case 4:
            return I74;
        case 5:
            return I75;
        default:
            throw std::runtime_error("Unsupported k for n=7.");
        }
    }

    default:
        throw std::runtime_error("Unsupported n.");
    }
}

inline double StableSqrtDifference(const double u1, const double u2, const double L1, const double L2) {
    return (u2 - u1) * (u2 + u1) / (L2 + L1);
}

inline double StableUOverPowerDifference(const double u1, const double u2, const double L1, const double L2,
                                         const int power, bool& reliable) {
    if (u1 == 0.0 || u2 == 0.0 || u1 * u2 < 0.0) {
        return u2 / std::pow(L2, power) - u1 / std::pow(L1, power);
    }

    const double logU = LogRatio(std::abs(u2), std::abs(u1));
    const double dL = StableSqrtDifference(u1, u2, L1, L2);
    const double logL = std::log1p(dL / L1);
    const double deltaLog = logU - static_cast<double>(power) * logL;
    const double logScale = std::abs(logU) + static_cast<double>(power) * std::abs(logL);
    if (std::abs(deltaLog) <= 1.0e-10 * logScale) {
        reliable = false;
    }

    const double value1 = u1 / std::pow(L1, power);
    return value1 * std::expm1(deltaLog);
}

inline double StableInversePowerDifference(const double u1, const double u2, const double L1, const double L2,
                                           const int power) {
    const double dL = StableSqrtDifference(u1, u2, L1, L2);
    const double logL = std::log1p(dL / L1);
    return std::pow(L1, -power) * std::expm1(-static_cast<double>(power) * logL);
}

inline double StableBaseMomentDiff(const double u1, const double u2, const double b, const int n, const int k,
                                   bool& reliable) {
    reliable = true;
    const double b2 = b * b;
    const double L1 = std::hypot(u1, b);
    const double L2 = std::hypot(u2, b);
    double D[4][6] = {};
    const bool even = (k % 2 == 0);

    if (even) {
        if (u1 > 0.0 && u2 > 0.0) {
            D[0][0] = LogRatio(u2 + L2, u1 + L1);
        } else if (u1 < 0.0 && u2 < 0.0) {
            D[0][0] = LogRatio(L1 - u1, L2 - u2);
        } else {
            D[0][0] = std::asinh(u2 / b) - std::asinh(u1 / b);
        }

        if (k >= 2) {
            double deltaUL;
            if (u1 * u2 > 0.0) {
                const double denominator = u2 * L2 + u1 * L1;
                if (denominator == 0.0) {
                    reliable = false;
                    deltaUL = u2 * L2 - u1 * L1;
                } else {
                    deltaUL = (u2 - u1) * (u2 + u1) * (u2 * u2 + u1 * u1 + b2) / denominator;
                }
            } else {
                deltaUL = u2 * L2 - u1 * L1;
            }
            D[0][2] = 0.5 * (deltaUL - b2 * D[0][0]);
        }
    } else {
        D[0][1] = StableSqrtDifference(u1, u2, L1, L2);
    }

    const int targetIndex = (n - 1) / 2;
    for (int index = 1; index <= targetIndex; ++index) {
        const int currentN = 2 * index + 1;
        const int power = currentN - 2;
        const int maxK = currentN == 3 ? 3 : (currentN == 5 ? 4 : 5);
        if (even) {
            const double inv = 1.0 / (static_cast<double>(power) * b2);
            bool endpointReliable = true;
            const double endpoint0 = StableUOverPowerDifference(u1, u2, L1, L2, power, endpointReliable);
            const double term0 = inv * endpoint0;
            const double term1 = static_cast<double>(currentN - 3) * inv * D[index - 1][0];
            D[index][0] = term0 + term1;
            reliable = reliable && endpointReliable;
            if (currentN > 3 && std::abs(D[index][0]) <= 1.0e-10 * (std::abs(term0) + std::abs(term1))) {
                reliable = false;
            }

            for (int currentK = 2; currentK <= std::min(k, maxK); currentK += 2) {
                const double first = D[index - 1][currentK - 2];
                const double second = b2 * D[index][currentK - 2];
                D[index][currentK] = std::fma(-b2, D[index][currentK - 2], first);
                if (std::abs(D[index][currentK]) <= 1.0e-10 * (std::abs(first) + std::abs(second))) {
                    reliable = false;
                }
            }
        } else {
            D[index][1] = -StableInversePowerDifference(u1, u2, L1, L2, power) / static_cast<double>(power);
            for (int currentK = 3; currentK <= std::min(k, maxK); currentK += 2) {
                const double first = D[index - 1][currentK - 2];
                const double second = b2 * D[index][currentK - 2];
                D[index][currentK] = std::fma(-b2, D[index][currentK - 2], first);
                if (std::abs(D[index][currentK]) <= 1.0e-10 * (std::abs(first) + std::abs(second))) {
                    reliable = false;
                }
            }
        }
    }

    const double result = D[targetIndex][k];
    reliable = reliable && std::isfinite(result);
    return result;
}

inline void StableBaseMomentDiffTable(const double u1, const double u2, const double b, const int n, const int maxK,
                                      double* output, bool& reliable) {
    reliable = true;
    const double b2 = b * b;
    const double L1 = std::hypot(u1, b);
    const double L2 = std::hypot(u2, b);
    double D[4][6] = {};

    if (u1 > 0.0 && u2 > 0.0) {
        D[0][0] = LogRatio(u2 + L2, u1 + L1);
    } else if (u1 < 0.0 && u2 < 0.0) {
        D[0][0] = LogRatio(L1 - u1, L2 - u2);
    } else {
        D[0][0] = std::asinh(u2 / b) - std::asinh(u1 / b);
    }

    if (maxK >= 1) {
        D[0][1] = StableSqrtDifference(u1, u2, L1, L2);
    }
    if (maxK >= 2) {
        double deltaUL;
        if (u1 * u2 > 0.0) {
            const double denominator = u2 * L2 + u1 * L1;
            if (denominator == 0.0) {
                reliable = false;
                deltaUL = u2 * L2 - u1 * L1;
            } else {
                deltaUL = (u2 - u1) * (u2 + u1) * (u2 * u2 + u1 * u1 + b2) / denominator;
            }
        } else {
            deltaUL = u2 * L2 - u1 * L1;
        }
        D[0][2] = 0.5 * (deltaUL - b2 * D[0][0]);
    }

    const int targetIndex = (n - 1) / 2;
    for (int index = 1; index <= targetIndex; ++index) {
        const int currentN = 2 * index + 1;
        const int power = currentN - 2;
        const double inv = 1.0 / (static_cast<double>(power) * b2);
        bool endpointReliable = true;
        const double endpoint0 = StableUOverPowerDifference(u1, u2, L1, L2, power, endpointReliable);
        const double term0 = inv * endpoint0;
        const double term1 = static_cast<double>(currentN - 3) * inv * D[index - 1][0];
        D[index][0] = term0 + term1;
        reliable = reliable && endpointReliable;
        if (currentN > 3 && std::abs(D[index][0]) <= 1.0e-10 * (std::abs(term0) + std::abs(term1))) {
            reliable = false;
        }

        if (maxK >= 1) {
            D[index][1] = -StableInversePowerDifference(u1, u2, L1, L2, power) / static_cast<double>(power);
        }

        const int supportedMaxK = currentN == 3 ? 3 : (currentN == 5 ? 4 : 5);
        for (int currentK = 2; currentK <= std::min(maxK, supportedMaxK); ++currentK) {
            const double first = D[index - 1][currentK - 2];
            const double second = b2 * D[index][currentK - 2];
            D[index][currentK] = std::fma(-b2, D[index][currentK - 2], first);
            if (std::abs(D[index][currentK]) <= 1.0e-10 * (std::abs(first) + std::abs(second))) {
                reliable = false;
            }
        }
    }

    for (int currentK = 0; currentK <= maxK; ++currentK) {
        output[currentK] = D[targetIndex][currentK];
        reliable = reliable && std::isfinite(output[currentK]);
    }
}

inline double BaseMomentGauss16(const double u1, const double u2, const double b, const int n, const int k) {
    static const double nodes[16] = {-0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                                     -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374,
                                     0.0950125098376374,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                                     0.7554044083550030,  0.8656312023878318,  0.9445750230732326,  0.9894009349916499};

    static const double weights[16] = {0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339,
                                       0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685,
                                       0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
                                       0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541};

    const double um = 0.5 * (u1 + u2);
    const double uh = 0.5 * (u2 - u1);

    double sum = 0.0;
    for (int i = 0; i < 16; ++i) {
        const double u = um + uh * nodes[i];
        const double f = std::pow(u, k) / std::pow(u * u + b * b, 0.5 * n);
        sum += weights[i] * f;
    }

    return uh * sum;
}

inline double BaseMomentDiff(const double u1, const double u2, const double bIn, const int n, const int k) {
    if (u1 == u2) {
        return 0.0;
    }

    double a = u1;
    double bnd = u2;
    double orderSign = 1.0;
    if (bnd < a) {
        std::swap(a, bnd);
        orderSign = -1.0;
    }

    double b = std::abs(bIn);
    if (b == 0.0) {
        b = std::numeric_limits<double>::min();
    }

    const bool sameSide = (a * bnd > 0.0);
    if (sameSide) {
        const double qmin = std::min(std::abs(a), std::abs(bnd));
        const double eta = b / std::max(qmin, std::numeric_limits<double>::min());
        if (eta < 1.0e-2) {
            return orderSign * BaseMomentSmallB(a, bnd, b, n, k);
        }
    }

    bool reliable = true;
    const double stable = StableBaseMomentDiff(a, bnd, b, n, k, reliable);
    if (reliable) {
        return orderSign * stable;
    }
    return orderSign * BaseMomentGauss16(a, bnd, b, n, k);
}

inline void BaseMomentDiffTable(const double u1, const double u2, const double bIn, const int n, const int maxK,
                                double* output) {
    if (u1 == u2) {
        std::fill(output, output + maxK + 1, 0.0);
        return;
    }

    double a = u1;
    double bnd = u2;
    double orderSign = 1.0;
    if (bnd < a) {
        std::swap(a, bnd);
        orderSign = -1.0;
    }

    double b = std::abs(bIn);
    if (b == 0.0) {
        b = std::numeric_limits<double>::min();
    }

    const bool sameSide = (a * bnd > 0.0);
    if (sameSide) {
        const double qmin = std::min(std::abs(a), std::abs(bnd));
        const double eta = b / std::max(qmin, std::numeric_limits<double>::min());
        if (eta < 1.0e-2) {
            for (int k = 0; k <= maxK; ++k) {
                output[k] = orderSign * BaseMomentSmallB(a, bnd, b, n, k);
            }
            return;
        }
    }

    bool reliable = true;
    StableBaseMomentDiffTable(a, bnd, b, n, maxK, output, reliable);
    if (reliable) {
        for (int k = 0; k <= maxK; ++k) {
            output[k] *= orderSign;
        }
        return;
    }

    for (int k = 0; k <= maxK; ++k) {
        output[k] = orderSign * BaseMomentGauss16(a, bnd, b, n, k);
    }
}

inline double RadialMoment(const double R2, const double R1, const double SinPsi, const double CosPsi, const int n,
                           const int p) {
    const double b = std::abs(SinPsi);
    const double c = ClampCos(CosPsi);
    const double u1 = R1 - c;
    const double u2 = R2 - c;

    double moments[6] = {};
    BaseMomentDiffTable(u1, u2, b, n, p, moments);
    double M = 0.0;
    for (int k = 0; k <= p; ++k) {
        const double coef = Binom(p, k) * std::pow(c, p - k);
        M += coef * moments[k];
    }
    return M;
}

inline double ZMoment(const double R2, const double R1, const double SinPsi, const double CosPsi, const int n,
                      const int p, const int q) {
    const double b = std::abs(SinPsi);
    const double c = ClampCos(CosPsi);
    const double u1 = R1 - c;
    const double u2 = R2 - c;

    double moments[6] = {};
    BaseMomentDiffTable(u1, u2, b, n, p + q, moments);
    double Z = 0.0;
    for (int k = 0; k <= p; ++k) {
        const double coefR = Binom(p, k) * std::pow(c, p - k);
        for (int m = 0; m <= q; ++m) {
            const double coefZ = Binom(q, m) * std::pow(c, m) * std::pow(-b * b, q - m);
            Z += coefR * coefZ * moments[k + m];
        }
    }
    return Z;
}

template <typename Computer>
inline Matrix14
EvaluateSurfaceKernel(const std::vector<double>& R2, const std::vector<double>& R1, const std::vector<double>& rho,
                      const Matrix14& FaiI, const Matrix14& LamdaI,
                      const double r, const double fai, const double lamda, const double Lons, const double Lone,
                      const double Lats, const double Late, const double dLon, const double dLat, const double Lon0s,
                      const double Lon0e, const double Lat0s, const double Lat0e, const Computer& compute) {
    constexpr int R = 14;
    constexpr int C = 14;

    Matrix14 Kernel;
    Matrix14 R2_interp, R1_interp, rho_interp;

    sphinterp::SphBilinearInterp(Lons, Lone, Lats, Late, dLon, dLat, R2, R1, rho, LamdaI, FaiI, Lon0s, Lon0e, Lat0s,
                                 Lat0e, R2_interp, R1_interp, rho_interp);

    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            const double weight = rho_interp(i, j);
            if (weight == 0.0 || R2_interp(i, j) == R1_interp(i, j)) {
                Kernel(i, j) = 0.0;
                continue;
            }

            double Tx, Ty, SinPsi, CosPsi;
            LocalGeometry(FaiI(i, j), LamdaI(i, j), fai, lamda, Tx, Ty, SinPsi, CosPsi);

            const double R2n = R2_interp(i, j) / r;
            const double R1n = R1_interp(i, j) / r;
            const double cosFai = std::cos(FaiI(i, j) * Phithreshold::Deg2Rad);

            Kernel(i, j) = weight * cosFai * compute(R2n, R1n, Tx, Ty, SinPsi, CosPsi);
        }
    }

    return Kernel;
}

inline double KV(double R2, double R1, double, double, double SinPsi, double CosPsi) {
    return RadialMoment(R2, R1, SinPsi, CosPsi, 1, 2);
}

inline double KVx(double R2, double R1, double Tx, double, double SinPsi, double CosPsi) {
    return Tx * RadialMoment(R2, R1, SinPsi, CosPsi, 3, 3);
}

inline double KVy(double R2, double R1, double, double Ty, double SinPsi, double CosPsi) {
    return Ty * RadialMoment(R2, R1, SinPsi, CosPsi, 3, 3);
}

inline double KVz(double R2, double R1, double, double, double SinPsi, double CosPsi) {
    return ZMoment(R2, R1, SinPsi, CosPsi, 3, 2, 1);
}

inline double KVxx(double R2, double R1, double Tx, double, double SinPsi, double CosPsi) {
    const double M32 = RadialMoment(R2, R1, SinPsi, CosPsi, 3, 2);
    const double M54 = RadialMoment(R2, R1, SinPsi, CosPsi, 5, 4);
    return 3.0 * Tx * Tx * M54 - M32;
}

inline double KVxy(double R2, double R1, double Tx, double Ty, double SinPsi, double CosPsi) {
    return 3.0 * Tx * Ty * RadialMoment(R2, R1, SinPsi, CosPsi, 5, 4);
}

inline double KVyy(double R2, double R1, double, double Ty, double SinPsi, double CosPsi) {
    const double M32 = RadialMoment(R2, R1, SinPsi, CosPsi, 3, 2);
    const double M54 = RadialMoment(R2, R1, SinPsi, CosPsi, 5, 4);
    return 3.0 * Ty * Ty * M54 - M32;
}

inline double KVzx(double R2, double R1, double Tx, double, double SinPsi, double CosPsi) {
    return 3.0 * Tx * ZMoment(R2, R1, SinPsi, CosPsi, 5, 3, 1);
}

inline double KVzy(double R2, double R1, double, double Ty, double SinPsi, double CosPsi) {
    return 3.0 * Ty * ZMoment(R2, R1, SinPsi, CosPsi, 5, 3, 1);
}

inline double KVzz(double R2, double R1, double, double, double SinPsi, double CosPsi) {
    const double M32 = RadialMoment(R2, R1, SinPsi, CosPsi, 3, 2);
    const double Z522 = ZMoment(R2, R1, SinPsi, CosPsi, 5, 2, 2);
    return 3.0 * Z522 - M32;
}

inline double KVxxx(double R2, double R1, double Tx, double, double SinPsi, double CosPsi) {
    const double M53 = RadialMoment(R2, R1, SinPsi, CosPsi, 5, 3);
    const double M75 = RadialMoment(R2, R1, SinPsi, CosPsi, 7, 5);
    return 15.0 * Tx * Tx * Tx * M75 - 9.0 * Tx * M53;
}

inline double KVxxy(double R2, double R1, double Tx, double Ty, double SinPsi, double CosPsi) {
    const double M53 = RadialMoment(R2, R1, SinPsi, CosPsi, 5, 3);
    const double M75 = RadialMoment(R2, R1, SinPsi, CosPsi, 7, 5);
    return 15.0 * Tx * Tx * Ty * M75 - 3.0 * Ty * M53;
}

inline double KVxxz(double R2, double R1, double Tx, double, double SinPsi, double CosPsi) {
    const double Z521 = ZMoment(R2, R1, SinPsi, CosPsi, 5, 2, 1);
    const double Z741 = ZMoment(R2, R1, SinPsi, CosPsi, 7, 4, 1);
    return 15.0 * Tx * Tx * Z741 - 3.0 * Z521;
}

inline double KVxyz(double R2, double R1, double Tx, double Ty, double SinPsi, double CosPsi) {
    return 15.0 * Tx * Ty * ZMoment(R2, R1, SinPsi, CosPsi, 7, 4, 1);
}

inline double KVyyx(double R2, double R1, double Tx, double Ty, double SinPsi, double CosPsi) {
    const double M53 = RadialMoment(R2, R1, SinPsi, CosPsi, 5, 3);
    const double M75 = RadialMoment(R2, R1, SinPsi, CosPsi, 7, 5);
    return 15.0 * Ty * Ty * Tx * M75 - 3.0 * Tx * M53;
}

inline double KVyyy(double R2, double R1, double, double Ty, double SinPsi, double CosPsi) {
    const double M53 = RadialMoment(R2, R1, SinPsi, CosPsi, 5, 3);
    const double M75 = RadialMoment(R2, R1, SinPsi, CosPsi, 7, 5);
    return 15.0 * Ty * Ty * Ty * M75 - 9.0 * Ty * M53;
}

inline double KVyyz(double R2, double R1, double, double Ty, double SinPsi, double CosPsi) {
    const double Z521 = ZMoment(R2, R1, SinPsi, CosPsi, 5, 2, 1);
    const double Z741 = ZMoment(R2, R1, SinPsi, CosPsi, 7, 4, 1);
    return 15.0 * Ty * Ty * Z741 - 3.0 * Z521;
}

inline double KVzzx(double R2, double R1, double Tx, double, double SinPsi, double CosPsi) {
    const double M53 = RadialMoment(R2, R1, SinPsi, CosPsi, 5, 3);
    const double Z732 = ZMoment(R2, R1, SinPsi, CosPsi, 7, 3, 2);
    return 15.0 * Tx * Z732 - 3.0 * Tx * M53;
}

inline double KVzzy(double R2, double R1, double, double Ty, double SinPsi, double CosPsi) {
    const double M53 = RadialMoment(R2, R1, SinPsi, CosPsi, 5, 3);
    const double Z732 = ZMoment(R2, R1, SinPsi, CosPsi, 7, 3, 2);
    return 15.0 * Ty * Z732 - 3.0 * Ty * M53;
}

inline double KVzzz(double R2, double R1, double, double, double SinPsi, double CosPsi) {
    const double Z521 = ZMoment(R2, R1, SinPsi, CosPsi, 5, 2, 1);
    const double Z723 = ZMoment(R2, R1, SinPsi, CosPsi, 7, 2, 3);
    return 15.0 * Z723 - 9.0 * Z521;
}
}

#define TESSEROID_SURFACE_KERNEL_FUNCTION(NAME, COMPUTER)                                                              \
    inline Matrix14 NAME(                                                                                              \
        const std::vector<double>& R2, const std::vector<double>& R1, const std::vector<double>& rho,                  \
        const Matrix14& FaiI, const Matrix14& LamdaI, const double r,                                                   \
        const double fai, const double lamda, const double Lons, const double Lone, const double Lats,                 \
        const double Late, const double dLon, const double dLat, const double Lon0s, const double Lon0e,               \
        const double Lat0s, const double Lat0e) {                                                                      \
        return surfacekernel_detail::EvaluateSurfaceKernel(R2, R1, rho, FaiI, LamdaI, r, fai, lamda, Lons, Lone, Lats, \
                                                           Late, dLon, dLat, Lon0s, Lon0e, Lat0s, Lat0e,               \
                                                           surfacekernel_detail::COMPUTER);                            \
    }

TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelV, KV)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVx, KVx)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVy, KVy)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVz, KVz)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVxx, KVxx)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVxy, KVxy)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVyy, KVyy)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVzx, KVzx)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVzy, KVzy)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVzz, KVzz)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVxxx, KVxxx)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVxxy, KVxxy)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVxxz, KVxxz)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVxyz, KVxyz)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVyyx, KVyyx)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVyyy, KVyyy)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVyyz, KVyyz)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVzzx, KVzzx)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVzzy, KVzzy)
TESSEROID_SURFACE_KERNEL_FUNCTION(Tesseroid_SurfaceKernelVzzz, KVzzz)

#undef TESSEROID_SURFACE_KERNEL_FUNCTION

}
