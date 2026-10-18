#include "mex.h"

#include "_TLI_estimate.h"

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::array<std::size_t, 20> gravity_output_order{{
    0, 1, 2, 3, 4, 5, 7, 6, 8, 9,
    10, 11, 12, 14, 13, 16, 17, 15, 18, 19
}};

struct MatrixView {
    const double* data = nullptr;
    mwSize rows = 0;
    mwSize columns = 0;

    double operator()(mwSize row, mwSize column) const {
        return data[row + rows * column];
    }
};

struct GatewayError : std::runtime_error {
    std::string id;
    GatewayError(const char* error_id, const std::string& message)
        : std::runtime_error(message), id(error_id) {}
};

void require(bool condition, const char* id, const std::string& message) {
    if (!condition) {
        throw GatewayError(id, message);
    }
}

MatrixView real_double_matrix(const mxArray* value, mwSize columns,
                              const char* argument_name) {
    require(mxIsDouble(value) && !mxIsComplex(value) && !mxIsSparse(value),
            "TesseroidLineIntegration:InvalidType",
            std::string(argument_name) + " must be a full, real double matrix.");
    require(mxGetNumberOfDimensions(value) == 2 && mxGetN(value) == columns,
            "TesseroidLineIntegration:InvalidSize",
            std::string(argument_name) + " must have exactly " +
                std::to_string(static_cast<unsigned long long>(columns)) + " columns.");
    return MatrixView{mxGetPr(value), mxGetM(value), mxGetN(value)};
}

double nonnegative_scalar(const mxArray* value, const char* argument_name) {
    require(mxIsDouble(value) && !mxIsComplex(value) &&
                mxGetNumberOfElements(value) == 1,
            "TesseroidLineIntegration:InvalidTolerance",
            std::string(argument_name) + " must be a real double scalar.");
    const double result = mxGetScalar(value);
    require(std::isfinite(result) && result >= 0.0,
            "TesseroidLineIntegration:InvalidTolerance",
            std::string(argument_name) + " must be finite and nonnegative.");
    return result;
}

std::string text_scalar(const mxArray* value) {
    require(mxIsChar(value),
            "TesseroidLineIntegration:InvalidMode",
            "The first input must be 'gravity' or 'magnetic'.");
    char* raw = mxArrayToUTF8String(value);
    require(raw != nullptr,
            "TesseroidLineIntegration:InvalidMode",
            "Unable to read the requested calculation mode.");
    std::string result(raw);
    mxFree(raw);
    return result;
}

bool logical_scalar(const mxArray* value) {
    require((mxIsLogical(value) || mxIsDouble(value)) &&
                mxGetNumberOfElements(value) == 1,
            "TesseroidLineIntegration:InvalidParallelFlag",
            "The parallel flag must be a logical or numeric scalar.");
    return mxIsLogical(value) ? mxIsLogicalScalarTrue(value)
                              : mxGetScalar(value) != 0.0;
}

Tesseroid::EstimateOptions make_options(double absolute_tolerance,
                                        double relative_tolerance) {
    Tesseroid::EstimateOptions options;
    options.integration.atol = absolute_tolerance;
    options.integration.rtol = relative_tolerance;
    options.gravitational_constant = 1.0;
    return options;
}

Tesseroid::Prism gravity_prism(const MatrixView& model, mwSize row) {
    constexpr double degree = Tesseroid::PI / 180.0;
    return Tesseroid::Prism{
        degree * model(row, 0), degree * model(row, 1),
        degree * model(row, 2), degree * model(row, 3),
        model(row, 4), model(row, 5), model(row, 6)
    };
}

Tesseroid::Point observation_point(const MatrixView& points, mwSize row) {
    constexpr double degree = Tesseroid::PI / 180.0;
    return Tesseroid::Point{
        degree * points(row, 0),
        Tesseroid::detail::normalize_longitude(degree * points(row, 1)),
        points(row, 2)
    };
}

struct Diagnostics {
    std::atomic<std::size_t> nonconverged{0};
    std::mutex mutex;
    std::string first_message;

    void inspect(const std::vector<Tesseroid::EstimateInfo>& result,
                 mwSize point_index, mwSize model_index) {
    for (const auto& item : result) {
        if (!std::isfinite(item.value) || !std::isfinite(item.error)) {
            std::ostringstream message;
            message << "Non-finite result for point "
                    << static_cast<unsigned long long>(point_index + 1)
                    << ", tesseroid "
                    << static_cast<unsigned long long>(model_index + 1)
                    << ", component " << Tesseroid::component_name(item.component)
                    << ". Estimated error=" << item.error
                    << ", value=" << item.value << ".";
            throw std::runtime_error(message.str());
        }
        if (!item.converged) {
            std::ostringstream message;
            message << "Point "
                    << static_cast<unsigned long long>(point_index + 1)
                    << ", tesseroid "
                    << static_cast<unsigned long long>(model_index + 1)
                    << ", component " << Tesseroid::component_name(item.component)
                    << ": estimated error=" << item.error
                    << ", value=" << item.value << ".";
            nonconverged.fetch_add(1, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(mutex);
            if (first_message.empty()) {
                first_message = message.str();
            }
        }
    }
    }

    void warn_if_needed() const {
        const std::size_t count = nonconverged.load(std::memory_order_relaxed);
        if (count != 0) {
            mexWarnMsgIdAndTxt(
                "TesseroidLineIntegration:Nonconvergence",
                "%llu component evaluation(s) did not meet the requested "
                "tolerance and were returned with finite values. First: %s",
                static_cast<unsigned long long>(count), first_message.c_str());
        }
    }
};

std::array<double, 20> gravity_at_point(
    const MatrixView& model, const MatrixView& points, mwSize point_index,
    const Tesseroid::EstimateOptions& options, Diagnostics& diagnostics) {
    const Tesseroid::Point point = observation_point(points, point_index);
    std::array<double, 20> sum{};
    std::array<double, 20> correction{};
    for (mwSize model_index = 0; model_index < model.rows; ++model_index) {
        const auto result = Tesseroid::estimate_all(
            gravity_prism(model, model_index), point, options);
        diagnostics.inspect(result, point_index, model_index);
        for (std::size_t component = 0; component < sum.size(); ++component) {
            Tesseroid::detail::compensated_add(
                result[component].value, sum[component], correction[component]);
        }
    }
    for (std::size_t component = 0; component < sum.size(); ++component) {
        sum[component] += correction[component];
    }
    return sum;
}

std::array<double, 10> magnetic_at_point(
    const MatrixView& model, const MatrixView& points, mwSize point_index,
    const Tesseroid::EstimateOptions& options) {
    const auto column = [&](mwSize index) {
        return model.data + model.rows * index;
    };
    const auto result = Tesseroid::work_detail::magnetic_coordinate(
        static_cast<int>(model.rows),
        column(0), column(1), column(2), column(3), column(4), column(5),
        column(6), column(7), column(8),
        points(point_index, 0), points(point_index, 1),
        points(point_index, 2), options);

    for (double value : result) {
        if (!std::isfinite(value)) {
            std::ostringstream message;
            message << "Non-finite magnetic result for point "
                    << static_cast<unsigned long long>(point_index + 1) << ".";
            throw std::runtime_error(message.str());
        }
    }
    return result;
}

template <std::size_t N, class Evaluator>
void evaluate_points(mwSize point_count, bool parallel, Evaluator&& evaluator,
                     const std::array<double*, N>& output) {
    std::atomic<bool> failed{false};
    std::mutex error_mutex;
    std::string error_message;

    auto evaluate_one = [&](mwSize point_index) {
        if (failed.load(std::memory_order_relaxed)) {
            return;
        }
        try {
            const auto value = evaluator(point_index);
            for (std::size_t component = 0; component < N; ++component) {
                output[component][point_index] = value[component];
            }
        } catch (const std::exception& error) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!failed.exchange(true)) {
                error_message = error.what();
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!failed.exchange(true)) {
                error_message = "Unknown C++ exception.";
            }
        }
    };

    if (parallel) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
        for (mwIndex point_index = 0;
             point_index < static_cast<mwIndex>(point_count); ++point_index) {
            evaluate_one(static_cast<mwSize>(point_index));
        }
#else
        throw std::runtime_error(
            "This MEX file was built without OpenMP support. Rebuild it with "
            "build_Tesseroid_LineIntegration_mex.");
#endif
    } else {
        for (mwSize point_index = 0; point_index < point_count; ++point_index) {
            evaluate_one(point_index);
            if (failed.load()) {
                break;
            }
        }
    }

    if (failed.load()) {
        throw std::runtime_error(error_message);
    }
}

template <std::size_t N>
std::array<double*, N> create_outputs(mxArray* output[], int output_count,
                                      mwSize point_count) {
    require(output_count == static_cast<int>(N),
            "TesseroidLineIntegration:OutputCount",
            "The MATLAB wrapper must request every documented output.");
    std::array<double*, N> pointers{};
    for (std::size_t component = 0; component < N; ++component) {
        output[component] = mxCreateDoubleMatrix(point_count, 1, mxREAL);
        pointers[component] = mxGetPr(output[component]);
    }
    return pointers;
}

} // namespace

MEXFUNCTION_LINKAGE void mexFunction(int output_count, mxArray* output[],
                                     int input_count, const mxArray* input[]) {
    try {
        require(input_count == 6,
                "TesseroidLineIntegration:InputCount",
                "Expected mode, TesseroidParameters, PointParameters, AT, RT, "
                "and a parallel flag.");
        const std::string mode = text_scalar(input[0]);
        const double absolute_tolerance = nonnegative_scalar(input[3], "AT");
        const double relative_tolerance = nonnegative_scalar(input[4], "RT");
        const bool parallel = logical_scalar(input[5]);
        const auto options = make_options(absolute_tolerance, relative_tolerance);

        if (mode == "gravity") {
            const MatrixView model =
                real_double_matrix(input[1], 7, "TesseroidParameters");
            const MatrixView points =
                real_double_matrix(input[2], 3, "PointParameters");
            const auto canonical_output =
                create_outputs<20>(output, output_count, points.rows);
            std::array<double*, 20> ordered_output{};
            Diagnostics diagnostics;
            for (std::size_t position = 0; position < ordered_output.size(); ++position) {
                ordered_output[gravity_output_order[position]] =
                    canonical_output[position];
            }
            evaluate_points<20>(
                points.rows, parallel,
                [&](mwSize point_index) {
                    return gravity_at_point(
                        model, points, point_index, options, diagnostics);
                },
                ordered_output);
            diagnostics.warn_if_needed();
            return;
        }

        if (mode == "magnetic") {
            const MatrixView model =
                real_double_matrix(input[1], 9, "TesseroidParameters");
            const MatrixView points =
                real_double_matrix(input[2], 3, "PointParameters");
            require(
                model.rows <=
                    static_cast<mwSize>(std::numeric_limits<int>::max()),
                "TesseroidLineIntegration:ModelTooLarge",
                "The number of magnetic tesseroids exceeds the C++ "
                "interface limit.");
            const auto pointers =
                create_outputs<10>(output, output_count, points.rows);
            evaluate_points<10>(
                points.rows, parallel,
                [&](mwSize point_index) {
                    return magnetic_at_point(
                        model, points, point_index, options);
                },
                pointers);
            return;
        }

        throw GatewayError("TesseroidLineIntegration:InvalidMode",
                           "Unknown mode '" + mode + "'. Use 'gravity' or 'magnetic'.");
    } catch (const GatewayError& error) {
        mexErrMsgIdAndTxt(error.id.c_str(), "%s", error.what());
    } catch (const std::exception& error) {
        mexErrMsgIdAndTxt("TesseroidLineIntegration:CalculationFailed",
                          "%s", error.what());
    } catch (...) {
        mexErrMsgIdAndTxt("TesseroidLineIntegration:CalculationFailed",
                          "Unknown native exception.");
    }
}
