// TLI (line-integration) MPI runtime benchmark for a single spherical prism.
//
// Model: the 1 km-above-surface case in "Validation by using Laplace's
// equation" (Test No. 8 in the manuscript)
//   latitude  [30, 50] deg
//   longitude [110, 130] deg
//   radius    [RE - 30 km, RE]
// Default observation radius: RE + 1 km (above the tesseroid).
//
// This is the TLI port of the original tesseroid RunTimeTest1 benchmark
// Differences from the original surface-quadrature version:
//   * Each component is evaluated once per point through
//     Tesseroid::estimate(...): the gauge-optimized line-integration kernel
//     with its own projection-split adaptive Gauss--Kronrod quadrature.
//   * The shell-minus-surroundings branch used by the original for points
//     inside the angular projection is not needed: the line-integration
//     representation treats interior, exterior, and polar-axis points
//     uniformly and carries the correct Poisson trace automatically.
//   * No Eigen dependency is required.
//
// The output CSV keeps the exact column layout of the original benchmark so
// that Summarize_RuntimeTest.m can be reused unchanged:
//   point_index, longitude_deg, latitude_deg, observation_height_m, radius_m,
//   inside_angular_bounds, distance_to_nearest_angular_boundary_deg,
//   <20 x component>_mean_s, all_components_mean_s, checksum
//
// Build example (TLI headers must be on the include path):
//   mpiicpx -O2 -std=c++17 -I. Test_RuntimeTest.cpp -o Test_RuntimeTest
//
// Run example:
//   mpirun -np 1800 ./Test_RuntimeTest
//
// NOTE: TLI evaluates every component with its own projection-split adaptive
// integration, so the all-components workload is about 20x one component.
// For a quick smoke test reduce the grid step
// or the timed repetitions below.

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "_TLI_estimate.h"
#include "_TLI_integralkernel.h"
#include "_integral.h"

namespace
{

// --------------------------- User parameters ---------------------------

constexpr double EARTH_RADIUS_M = 6371200.0;
// Height relative to EARTH_RADIUS_M. It may be positive, zero, or negative.
// Manuscript cases:
//   Test No. 8:  +1000 m
//   Test No. 9: -31000 m
//   Test No. 10: -15000 m
// Default runtime example:
//   1 km above the outer surface: +1000 m
constexpr double OBSERVATION_HEIGHT_M = 1000.0;

constexpr double PRISM_LAT_MIN_DEG = 30.0;
constexpr double PRISM_LAT_MAX_DEG = 50.0;
constexpr double PRISM_LON_MIN_DEG = 110.0;
constexpr double PRISM_LON_MAX_DEG = 130.0;
constexpr double PRISM_R_MIN_M = EARTH_RADIUS_M - 30000.0;
constexpr double PRISM_R_MAX_M = EARTH_RADIUS_M;
constexpr double PRISM_DENSITY = 2670.0;

// Runtime-test grid over the manuscript's Laplace/Poisson validation domain.
// A 0.1-degree spacing gives 401 x 401 = 160801 points.
constexpr double GRID_LAT_MIN_DEG = 20.05;
constexpr double GRID_LAT_MAX_DEG = 60.05;
constexpr double GRID_LAT_STEP_DEG = 0.1;
constexpr double GRID_LON_MIN_DEG = 100.05;
constexpr double GRID_LON_MAX_DEG = 140.05;
constexpr double GRID_LON_STEP_DEG = 0.1;

constexpr int WARMUP_REPETITIONS = 2;
constexpr int TIMED_REPETITIONS = 10;

constexpr double ABS_TOL = 1.0e-16;
constexpr double REL_TOL = 1.0e-16;

// The runtime is geometry/kernel-dominated; the value of G only scales the
// checksum, not the measured time.  TLI is length-unit invariant.
constexpr double MODEL_GRAVITATIONAL_CONSTANT = 1.0;

constexpr const char* OUTPUT_FILE =
    "RuntimeTest_LaplacePoisson_HeightPlus1km.csv";

// ----------------------------------------------------------------------

constexpr int COMPONENT_COUNT = 20;

// Work-compatible component labels in the original Result_TFG.dat order.
// TLI's component(name) maps the symmetric aliases (Vzx->Vxz, Vzy->Vyz,
// Vyyx->Vxyy, Vzzx->Vxzz, Vzzy->Vyzz) to the canonical tensor order.
const std::array<const char*, COMPONENT_COUNT> COMPONENT_NAMES = {
    "V",
    "Vx", "Vy", "Vz",
    "Vxx", "Vxy", "Vyy", "Vzx", "Vzy", "Vzz",
    "Vxxx", "Vxxy", "Vxxz", "Vxyz", "Vyyx",
    "Vyyy", "Vyyz", "Vzzx", "Vzzy", "Vzzz"
};

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double DEG_TO_RAD = PI / 180.0;

struct GridPoint
{
    double latitude = 0.0;
    double longitude = 0.0;
};

struct RuntimeRecord
{
    double latitude = 0.0;
    double longitude = 0.0;
    std::array<double, COMPONENT_COUNT> meanSeconds{};
    double allComponentsMeanSeconds = 0.0;
    double checksum = 0.0;
};

void synchronizedProgress(int rank, const std::string& message)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "[all ranks] " << message << '\n' << std::flush;
    }
}

void mpiDistribute(int total, std::vector<int>& counts,
                   std::vector<int>& displacements)
{
    const int size = static_cast<int>(counts.size());
    for (int p = 0; p < size; ++p) {
        counts[p] = total / size + (p < total % size ? 1 : 0);
        displacements[p] = p == 0 ? 0 : displacements[p - 1] + counts[p - 1];
    }
}

std::vector<double> makeAxis(double begin, double end, double step)
{
    if (!(step > 0.0) || end < begin) {
        throw std::runtime_error("Invalid grid-axis parameters.");
    }

    const long long count =
        static_cast<long long>(std::floor((end - begin) / step + 0.5)) + 1;

    std::vector<double> axis;
    axis.reserve(static_cast<std::size_t>(count));
    for (long long i = 0; i < count; ++i) {
        axis.push_back(begin + static_cast<double>(i) * step);
    }
    return axis;
}

// Evaluate one gravitational component with the TLI line-integration kernel.
double evaluateOneComponent(
    const char* componentName,
    double latitude,
    double longitude,
    double radius)
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
    options.integration.atol = ABS_TOL;
    options.integration.rtol = REL_TOL;
    options.gravitational_constant = MODEL_GRAVITATIONAL_CONSTANT;

    EstimateInfo info;
    const double value = estimate(component(componentName), prism, point,
                                  options, &info);

    if (!std::isfinite(value)) {
        throw std::runtime_error(
            std::string("TLI estimate returned a non-finite value for ") +
            componentName + " at (" +
            std::to_string(latitude) + ", " + std::to_string(longitude) + ")");
    }

    return value;
}

RuntimeRecord benchmarkPoint(const GridPoint& point)
{
    const double radius = EARTH_RADIUS_M + OBSERVATION_HEIGHT_M;

    RuntimeRecord record;
    record.latitude = point.latitude;
    record.longitude = point.longitude;

    volatile double sink = 0.0;

    for (int component = 0; component < COMPONENT_COUNT; ++component) {
        for (int repeat = 0; repeat < WARMUP_REPETITIONS; ++repeat) {
            sink += evaluateOneComponent(
                COMPONENT_NAMES[component],
                point.latitude,
                point.longitude,
                radius);
        }

        const double start = MPI_Wtime();
        for (int repeat = 0; repeat < TIMED_REPETITIONS; ++repeat) {
            sink += evaluateOneComponent(
                COMPONENT_NAMES[component],
                point.latitude,
                point.longitude,
                radius);
        }
        const double elapsed = MPI_Wtime() - start;
        record.meanSeconds[component] =
            elapsed / static_cast<double>(TIMED_REPETITIONS);
    }

    record.allComponentsMeanSeconds = 0.0;
    for (double seconds : record.meanSeconds) {
        record.allComponentsMeanSeconds += seconds;
    }

    record.checksum = sink;
    return record;
}

void writeCsv(
    const std::vector<GridPoint>& points,
    const std::vector<double>& gathered,
    int valuesPerPoint)
{
    std::ofstream output(OUTPUT_FILE);
    if (!output) {
        throw std::runtime_error(
            std::string("Cannot create output file: ") + OUTPUT_FILE);
    }

    output << std::setprecision(17);
    output << "point_index,longitude_deg,latitude_deg,"
           << "observation_height_m,radius_m,"
           << "inside_angular_bounds,"
           << "distance_to_nearest_angular_boundary_deg";
    for (const char* name : COMPONENT_NAMES) {
        output << ',' << name << "_mean_s";
    }
    output << ",all_components_mean_s,checksum\n";

    const double radius = EARTH_RADIUS_M + OBSERVATION_HEIGHT_M;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const std::size_t offset = i * static_cast<std::size_t>(valuesPerPoint);
        const bool insideAngularBounds =
            points[i].latitude > PRISM_LAT_MIN_DEG &&
            points[i].latitude < PRISM_LAT_MAX_DEG &&
            points[i].longitude > PRISM_LON_MIN_DEG &&
            points[i].longitude < PRISM_LON_MAX_DEG;

        const double boundaryDistance = std::min(
            std::min(
                std::abs(points[i].latitude - PRISM_LAT_MIN_DEG),
                std::abs(points[i].latitude - PRISM_LAT_MAX_DEG)),
            std::min(
                std::abs(points[i].longitude - PRISM_LON_MIN_DEG),
                std::abs(points[i].longitude - PRISM_LON_MAX_DEG)));

        output << i << ','
               << points[i].longitude << ','
               << points[i].latitude << ','
               << OBSERVATION_HEIGHT_M << ','
               << radius << ','
               << (insideAngularBounds ? 1 : 0) << ','
               << boundaryDistance;

        for (int component = 0; component < COMPONENT_COUNT; ++component) {
            output << ',' << gathered[offset + component];
        }
        output << ',' << gathered[offset + COMPONENT_COUNT]
               << ',' << gathered[offset + COMPONENT_COUNT + 1]
               << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int processCount = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    int exitCode = 0;

    try {
        const std::vector<double> latitudes =
            makeAxis(GRID_LAT_MIN_DEG, GRID_LAT_MAX_DEG, GRID_LAT_STEP_DEG);
        const std::vector<double> longitudes =
            makeAxis(GRID_LON_MIN_DEG, GRID_LON_MAX_DEG, GRID_LON_STEP_DEG);

        std::vector<GridPoint> points;
        points.reserve(latitudes.size() * longitudes.size());
        for (double latitude : latitudes) {
            for (double longitude : longitudes) {
                points.push_back(GridPoint{latitude, longitude});
            }
        }

        std::vector<int> counts(processCount, 0);
        std::vector<int> displacements(processCount, 0);
        mpiDistribute(static_cast<int>(points.size()), counts, displacements);

        const int localCount = counts[rank];
        const int localBegin = displacements[rank];
        constexpr int valuesPerPoint = COMPONENT_COUNT + 2;

        std::vector<double> localResults(
            static_cast<std::size_t>(localCount) * valuesPerPoint, 0.0);

        if (rank == 0) {
            std::cout
                << "TLI line-integration runtime benchmark\n"
                << "Laplace/Poisson validation case 1 "
                   "(manuscript Test No. 8, height +1 km)\n"
                << "MPI processes: " << processCount << '\n'
                << "Grid points: " << points.size()
                << " (" << longitudes.size() << " longitude x "
                << latitudes.size() << " latitude)\n"
                << "Warm-up repetitions: " << WARMUP_REPETITIONS << '\n'
                << "Timed repetitions: " << TIMED_REPETITIONS << '\n'
                << "Observation height relative to RE: "
                << OBSERVATION_HEIGHT_M << " m\n"
                << "Observation radius: "
                << (EARTH_RADIUS_M + OBSERVATION_HEIGHT_M) << " m\n"
                << "Output: " << OUTPUT_FILE << "\n\n";
        }

        synchronizedProgress(rank, "initial synchronization complete; starting local benchmark");
        const double benchmarkStart = MPI_Wtime();

        for (int localIndex = 0; localIndex < localCount; ++localIndex) {
            const int globalIndex = localBegin + localIndex;
            const RuntimeRecord record = benchmarkPoint(points[globalIndex]);
            const std::size_t offset =
                static_cast<std::size_t>(localIndex) * valuesPerPoint;

            for (int component = 0; component < COMPONENT_COUNT; ++component) {
                localResults[offset + component] =
                    record.meanSeconds[component];
            }
            localResults[offset + COMPONENT_COUNT] =
                record.allComponentsMeanSeconds;
            localResults[offset + COMPONENT_COUNT + 1] = record.checksum;
        }

        const double localElapsed = MPI_Wtime() - benchmarkStart;

        std::vector<int> resultCounts(processCount, 0);
        std::vector<int> resultDisplacements(processCount, 0);
        for (int source = 0; source < processCount; ++source) {
            resultCounts[source] = counts[source] * valuesPerPoint;
            resultDisplacements[source] =
                displacements[source] * valuesPerPoint;
        }

        synchronizedProgress(rank, "local benchmark complete; entering MPI_Gatherv");

        std::vector<double> gatheredResults;
        if (rank == 0) {
            gatheredResults.resize(points.size() * valuesPerPoint);
        }

        MPI_Gatherv(
            localResults.data(),
            resultCounts[rank],
            MPI_DOUBLE,
            rank == 0 ? gatheredResults.data() : nullptr,
            resultCounts.data(),
            resultDisplacements.data(),
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD);

        double maximumElapsed = 0.0;
        MPI_Reduce(
            &localElapsed,
            &maximumElapsed,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0,
            MPI_COMM_WORLD);

        synchronizedProgress(rank, "MPI_Gatherv and MPI_Reduce complete");

        if (rank == 0) {
            std::cout << "[rank 0] writing CSV\n" << std::flush;
            writeCsv(points, gatheredResults, valuesPerPoint);
            std::cout << "\nBenchmark completed in "
                      << maximumElapsed
                      << " s (maximum rank computation time).\n"
                      << "CSV written to " << OUTPUT_FILE << '\n'
                      << std::flush;
        }

        synchronizedProgress(rank, "CSV output complete; entering MPI_Finalize");
    }
    catch (const std::exception& error) {
        std::cerr << "Rank " << rank << " error: " << error.what() << '\n';
        exitCode = 1;
    }

    MPI_Finalize();

    if (rank == 0) {
        std::cout << "finished MPI_Finalize\n" << std::flush;
    }
    return exitCode;
}
