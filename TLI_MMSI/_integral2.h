#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <vector>

namespace Tesseroid {

template <std::size_t R, std::size_t C> struct FixedMatrix {
    std::array<double, R * C> values{};

    double& operator()(std::size_t row, std::size_t col) {
        return values[row * C + col];
    }

    const double& operator()(std::size_t row, std::size_t col) const {
        return values[row * C + col];
    }

    void fill(double value) { values.fill(value); }
};

using Matrix14 = FixedMatrix<14, 14>;

struct Gauss3Kronrod7Const {
    const std::array<double, 7> low_weights{{
        0.0, 5.0 / 9.0, 0.0, 8.0 / 9.0, 0.0, 5.0 / 9.0, 0.0}};
    const std::array<double, 7> high_weights{{
        0.1046562260264672, 0.2684880898683334,
        0.4013974147759622, 0.4509165386584744,
        0.4013974147759622, 0.2684880898683334,
        0.1046562260264672}};
    const std::array<double, 14> nodes{{
        0.009877182822994962, 0.05635083268962915,
        0.14143906266329936, 0.25, 0.35856093733670064,
        0.44364916731037085, 0.490122817177005,
        0.509877182822995, 0.5563508326896291,
        0.6414390626632993, 0.75, 0.8585609373367007,
        0.9436491673103709, 0.990122817177005}};
};

struct SubRectangle {
    double q = 0.0;
    double e = 0.0;
    double L = 0.0;
    double R = 0.0;
    double B = 0.0;
    double T = 0.0;
    double adjer = 0.0;
};

struct Integral2Info {
    double errbnd = 0.0;
    bool minRectWarn = false;
    bool maxNFEWarn = false;
    int NFE = 0;
};

struct TensorEvalResult {
    bool subdivided = true;
    std::array<double, 4> qsub{};
    std::array<double, 4> esub{};
};

inline double sum4(const std::array<double, 4>& values) {
    return values[0] + values[1] + values[2] + values[3];
}

class RectangleQueue {
  private:
    struct Compare {
        bool operator()(const SubRectangle& a, const SubRectangle& b) const {
            return a.adjer < b.adjer;
        }
    };

    std::priority_queue<SubRectangle, std::vector<SubRectangle>, Compare> queue_;
    double err_ok_ = 0.0;
    double sum_adjer_ = 0.0;

    void pushRect(double q, double e, double L, double R, double B,
                  double T, double adjer) {
        queue_.push(SubRectangle{q, e, L, R, B, T, adjer});
        sum_adjer_ += adjer;
    }

  public:
    bool empty() const { return queue_.empty(); }
    double errBnd() const { return err_ok_ + sum_adjer_; }
    void addErrOk(double value) { err_ok_ += value; }

    SubRectangle pop() {
        const SubRectangle rect = queue_.top();
        queue_.pop();
        sum_adjer_ -= rect.adjer;
        return rect;
    }

    void saveRectInfo(const std::array<double, 4>& qsub,
                      const std::array<double, 4>& esub,
                      double thetaL, double thetaR, double phiB, double phiT,
                      double tolerance, double area, double eps100,
                      double adjust) {
        const double thetaM = 0.5 * (thetaL + thetaR);
        const double phiM = 0.5 * (phiB + phiT);
        const double subArea = 0.25 * (thetaR - thetaL) * (phiT - phiB);
        double localTolerance = tolerance * subArea / area;
        localTolerance = std::max(
            std::abs(localTolerance), eps100 * std::abs(sum4(qsub)));

        auto save = [&](std::size_t i, double L, double R, double B, double T) {
            const double adjustedError = adjust * esub[i];
            if (adjustedError > localTolerance)
                pushRect(qsub[i], esub[i], L, R, B, T, adjustedError);
            else
                err_ok_ += adjustedError;
        };

        save(0, thetaL, thetaM, phiB, phiM);
        save(1, thetaM, thetaR, phiB, phiM);
        save(2, thetaL, thetaM, phiM, phiT);
        save(3, thetaM, thetaR, phiM, phiT);
    }
};

template <class Func>
inline TensorEvalResult PartitionIntegralQuad(
    double thetaL, double thetaR, double phiB, double phiT,
    bool& firstFunEval, int& nfe, double xmin, double xmax,
    double ymin, double ymax, const Func& function) {
    TensorEvalResult result;
    constexpr double tiny = 100.0 * std::numeric_limits<double>::epsilon();
    if (std::abs(thetaR - thetaL) <= tiny ||
        std::abs(phiT - phiB) <= tiny || std::abs(xmax - xmin) <= tiny ||
        std::abs(ymax - ymin) <= tiny) {
        result.subdivided = false;
        return result;
    }

    const Gauss3Kronrod7Const rule{};
    const double dtheta = thetaR - thetaL;
    const double dphi = phiT - phiB;
    const double dydt = ymax - ymin;
    std::array<double, 14> theta{}, phi{}, x{}, y{};
    Matrix14 X, Y;

    for (std::size_t i = 0; i < 14; ++i) {
        theta[i] = thetaL + rule.nodes[i] * dtheta;
        phi[i] = phiB + rule.nodes[i] * dphi;
        x[i] = 0.5 * (xmax + xmin) +
               0.5 * (xmax - xmin) * std::cos(theta[i]);
        y[i] = ymin + (0.5 + 0.5 * std::cos(phi[i])) * dydt;
    }
    for (std::size_t i = 0; i < 14; ++i) {
        for (std::size_t j = 0; j < 14; ++j) {
            X(i, j) = x[j];
            Y(i, j) = y[i];
        }
    }

    Matrix14 values = function(X, Y);
    ++nfe;
    firstFunEval = false;
    const double scale = 0.0625 * dtheta * dphi;

    for (std::size_t phiHalf = 0; phiHalf < 2; ++phiHalf) {
        for (std::size_t thetaHalf = 0; thetaHalf < 2; ++thetaHalf) {
            double high = 0.0;
            double low = 0.0;
            for (std::size_t i = 0; i < 7; ++i) {
                const std::size_t row = phiHalf * 7 + i;
                const double phiJacobian = dydt * std::sin(phi[row]);
                for (std::size_t j = 0; j < 7; ++j) {
                    const std::size_t col = thetaHalf * 7 + j;
                    const double jacobian = 0.25 * (xmax - xmin) *
                                            phiJacobian * std::sin(theta[col]);
                    const double weightedValue = values(row, col) * jacobian;
                    high += rule.high_weights[i] * rule.high_weights[j] *
                            weightedValue;
                    low += rule.low_weights[i] * rule.low_weights[j] *
                           weightedValue;
                }
            }
            const std::size_t quadrant = phiHalf * 2 + thetaHalf;
            result.qsub[quadrant] = high * scale;
            result.esub[quadrant] = std::abs(low * scale -
                                                    result.qsub[quadrant]);
        }
    }
    return result;
}

template <class Func>
inline double integral2t(const Func& function, double xmin, double xmax,
                         double ymin, double ymax, double absoluteTolerance = 1e-10,
                         double relativeTolerance = 1e-6, int maxFunEvals = 1000,
                         Integral2Info* info = nullptr) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    bool firstFunEval = true;
    int nfe = 0;
    const double thetaL = 0.0;
    const double thetaR = pi;
    const double phiB = 0.0;
    const double phiT = pi;
    const double area = pi * pi;

    auto first = PartitionIntegralQuad(
        thetaL, thetaR, phiB, phiT, firstFunEval, nfe,
        xmin, xmax, ymin, ymax, function);
    double integral = sum4(first.qsub);
    const double eps100 = 100.0 * std::numeric_limits<double>::epsilon();
    relativeTolerance = std::max(relativeTolerance, eps100);
    const double relativeEighth = std::max(relativeTolerance / 8.0, eps100);
    const double absoluteEighth = absoluteTolerance / 8.0;
    double tolerance = std::max(absoluteEighth,
                                relativeEighth * std::abs(integral));
    double adjust = 1.0;

    RectangleQueue manager;
    manager.saveRectInfo(first.qsub, first.esub, thetaL, thetaR, phiB, phiT,
                         tolerance, area, eps100, adjust);
    bool minRectWarn = false;
    bool maxNFEWarn = false;

    while (!manager.empty() && manager.errBnd() > tolerance) {
        const SubRectangle rect = manager.pop();
        auto next = PartitionIntegralQuad(
            rect.L, rect.R, rect.B, rect.T, firstFunEval, nfe,
            xmin, xmax, ymin, ymax, function);
        if (!next.subdivided) {
            minRectWarn = true;
            manager.addErrOk(rect.adjer);
            continue;
        }

        const double newIntegral = sum4(next.qsub);
        adjust = rect.e > 0.0
                     ? std::min(1.0, std::abs(rect.q - newIntegral) / rect.e)
                     : 1.0;
        integral += newIntegral - rect.q;
        tolerance = std::max(absoluteEighth,
                             relativeEighth * std::abs(integral));
        manager.saveRectInfo(next.qsub, next.esub, rect.L, rect.R,
                             rect.B, rect.T, tolerance, area, eps100, adjust);
        if (nfe >= maxFunEvals) {
            maxNFEWarn = true;
            break;
        }
    }

    if (info != nullptr) {
        info->errbnd = manager.errBnd();
        info->minRectWarn = minRectWarn;
        info->maxNFEWarn = maxNFEWarn;
        info->NFE = nfe;
    }
    return integral;
}

}
