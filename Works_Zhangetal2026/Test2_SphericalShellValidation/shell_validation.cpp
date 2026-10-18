// TLI-MMSI spherical-shell validation.
// Reproduces Tests No. 5-7 of Zhang et al. (2025), "Accurate and Stable
// Computation of Gravitational Potential and Its Derivatives up to Third
// Order for Spherical Prisms".
//
// A complete 30-km-thick spherical shell (r1=6341.2 km, r2=6371.2 km,
// G=1, rho=1), discretized as a 1 deg x 1 deg cell-centred TLI-MMSI grid,
// is evaluated on the 1-deg global observation grid (64,800 points) at
//   Test 5: r = R_E + 1 km   (outside)
//   Test 6: r = R_E - 31 km  (inside)
//   Test 7: r = R_E - 15 km  (within the shell)
// and compared against the closed-form complete-spherical-shell expressions
// of the manuscript's Appendix B.
//
// Outputs:
//   results/shell_summary.csv                       error statistics
//   results/ResultsPan_10x10Split_1km.dat           Test 5 point data
//   results/ResultsPan_10x10Split_-31km.dat         Test 6 point data
//   results/ResultsPan_10x10Split_-15km.dat         Test 7 point data
// The .dat files use the manuscript Work binary layout:
//   int32 NumCoordinates, then 23 double arrays of length NumCoordinates:
//   Longitude, Latitude, Radius, V, Vx, Vy, Vz, Vxx, Vxy, Vyy, Vzx, Vzy,
//   Vzz, Vxxx, Vxxy, Vxxz, Vxyz, Vyyx, Vyyy, Vyyz, Vzzx, Vzzy, Vzzz.
//
// Build (Intel oneAPI, Windows):
//   call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
//   icpx -O2 -std=c++17 -DTLI_ENABLE_MPI
//        -I..\..\TLI_MMSI -I..\..\C++Codes
//        -I"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\include"
//        shell_validation.cpp -o shell_validation.exe
//        -L"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\lib" -limpi -limpicxx
// Run:
//   mpiexec -n 8 shell_validation.exe [stride]

#include "_tesseroid_MSI.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
inline int _mkdir(const char* path) { return mkdir(path, 0755); }
#endif
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr double R_E = 6371.2;          // reference Earth radius (km)
constexpr double SHELL_R1 = R_E - 30.0; // 6341.2 km
constexpr double SHELL_R2 = R_E;        // 6371.2 km

// Cell-centred 1 deg x 1 deg global grid: 360 lon x 180 lat = 64,800 cells.
constexpr double MODEL_LONS = -179.5;
constexpr double MODEL_LONE = 179.5;
constexpr double MODEL_LATS = -89.5;
constexpr double MODEL_LATE = 89.5;
constexpr double MODEL_DLON = 1.0;
constexpr double MODEL_DLAT = 1.0;
constexpr int MODEL_NLON = 360;
constexpr int MODEL_NLAT = 180;

// 1-deg observation grid: 360 lon x 180 lat = 64,800 points.
constexpr double OBS_LON_FIRST = -179.5;
constexpr double OBS_LAT_FIRST = -89.5;
constexpr int OBS_NLON = 360;
constexpr int OBS_NLAT = 180;
constexpr int OBS_COUNT = OBS_NLON * OBS_NLAT;

constexpr double PI = 3.141592653589793238462643383279502884;

// Canonical TLI 20-component order:
// 0 V, 1 Vx, 2 Vy, 3 Vz,
// 4 Vxx, 5 Vxy, 6 Vxz, 7 Vyy, 8 Vyz, 9 Vzz,
// 10 Vxxx, 11 Vxxy, 12 Vxxz, 13 Vxyy, 14 Vxyz, 15 Vxzz,
// 16 Vyyy, 17 Vyyz, 18 Vyzz, 19 Vzzz.
const char* component_name(int k) {
    static const char* names[] = {
        "V", "Vx", "Vy", "Vz", "Vxx", "Vxy", "Vxz", "Vyy", "Vyz", "Vzz",
        "Vxxx", "Vxxy", "Vxxz", "Vxyy", "Vxyz", "Vxzz", "Vyyy", "Vyyz",
        "Vyzz", "Vzzz"};
    return names[k];
}

int component_order(int k) {
    return k == 0 ? 0 : (k <= 3 ? 1 : (k <= 9 ? 2 : 3));
}

// Canonical TLI index -> manuscript Work order index (same permutation as the
// Work-file store() helper): V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz, Vxxx,Vxxy,
// Vxxz,Vxyz,Vyyx,Vyyy,Vyyz,Vzzx,Vzzy,Vzzz.
constexpr int work_order[20] = {
    0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10, 11, 12, 14, 13, 16, 17, 15, 18, 19};

// Closed-form complete-spherical-shell reference (Appendix B), km, G=1, rho=1.
struct ShellReference {
    double V, Vz, Vxx, Vzz, Vxxz, Vzzz;
};

inline ShellReference shell_reference(double r) {
    const double d2 = (SHELL_R2 - SHELL_R1) * (SHELL_R2 + SHELL_R1);
    const double d3 = (SHELL_R2 - SHELL_R1) *
                      (SHELL_R2 * SHELL_R2 + SHELL_R2 * SHELL_R1 +
                       SHELL_R1 * SHELL_R1);
    if (r >= SHELL_R2) {
        return {4.0 * PI * d3 / (3.0 * r),
                -4.0 * PI * d3 / (3.0 * r * r),
                -4.0 * PI * d3 / (3.0 * r * r * r),
                8.0 * PI * d3 / (3.0 * r * r * r),
                4.0 * PI * d3 / (r * r * r * r),
                -8.0 * PI * d3 / (r * r * r * r)};
    }
    if (r <= SHELL_R1) {
        return {2.0 * PI * d2, 0.0, 0.0, 0.0, 0.0, 0.0};
    }
    const double dr = (r - SHELL_R1) *
                      (r * r + r * SHELL_R1 + SHELL_R1 * SHELL_R1);
    return {2.0 * PI * (d2 - (r - SHELL_R1) * (r - SHELL_R1) *
                                (r + 2.0 * SHELL_R1) / (3.0 * r)),
            -4.0 * PI * dr / (3.0 * r * r),
            -4.0 * PI * dr / (3.0 * r * r * r),
            -4.0 * PI / 3.0 * (1.0 + 2.0 * SHELL_R1 * SHELL_R1 * SHELL_R1 /
                                          (r * r * r)),
            -4.0 * PI * SHELL_R1 * SHELL_R1 * SHELL_R1 / (r * r * r * r),
            8.0 * PI * SHELL_R1 * SHELL_R1 * SHELL_R1 / (r * r * r * r)};
}

inline std::array<double, 20> analytical(const ShellReference& s) {
    std::array<double, 20> a{};
    a[0] = s.V;
    a[3] = s.Vz;
    a[4] = s.Vxx;   // Vxx
    a[7] = s.Vxx;   // Vyy = Vxx
    a[9] = s.Vzz;
    a[12] = s.Vxxz;
    a[17] = s.Vxxz; // Vyyz = Vxxz
    a[19] = s.Vzzz;
    return a;
}

struct Stats {
    double min = std::numeric_limits<double>::infinity();
    double max = 0.0;
    double sum = 0.0;
    double sumsq = 0.0;
    std::size_t count = 0;
    std::size_t nan_count = 0;
};

inline void accumulate(Stats& s, double numerical, double reference) {
    const double err = std::abs(numerical - reference);
    if (!std::isfinite(err)) {
        ++s.nan_count;
        return;
    }
    s.min = std::min(s.min, err);
    s.max = std::max(s.max, err);
    s.sum += err;
    s.sumsq += err * err;
    ++s.count;
}

inline int local_count_int_for_rank(int rank, int total, int processes) {
    return total / processes + (rank < total % processes ? 1 : 0);
}

} // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int processes = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    int stride = 1;
    if (argc > 1)
        stride = std::max(1, std::atoi(argv[1]));

    int first_test = 0;
    int last_test = 2;
    if (argc > 2) {
        const int requested_test = std::atoi(argv[2]);
        if (requested_test < 5 || requested_test > 7) {
            if (rank == 0)
                std::fprintf(stderr,
                             "Usage: shell_validation.exe [stride] [test 5|6|7]\n");
            MPI_Finalize();
            return 2;
        }
        first_test = last_test = requested_test - 5;
    }
    if (argc > 3) {
        if (rank == 0)
            std::fprintf(stderr,
                         "Usage: shell_validation.exe [stride] [test 5|6|7]\n");
        MPI_Finalize();
        return 2;
    }

    // The executable is also intended to be runnable directly, without the
    // run.bat wrapper.  Create the output directory before any rank writes a
    // slice; otherwise a direct MPI launch silently loses all point data.
    if (rank == 0 && _mkdir("results") != 0 && errno != EEXIST) {
        std::fprintf(stderr, "Error: cannot create results directory\n");
        MPI_Abort(MPI_COMM_WORLD, 3);
    }
    MPI_Barrier(MPI_COMM_WORLD);

// ---------- Build the 1 deg x 1 deg cell-centred shell model ----------
    const std::size_t model_count =
        static_cast<std::size_t>(MODEL_NLON) * static_cast<std::size_t>(MODEL_NLAT);
    std::vector<double> model_r2(model_count, SHELL_R2);
    std::vector<double> model_r1(model_count, SHELL_R1);
    std::vector<double> model_rho(model_count, 1.0);

    const Tesseroid::tli_mmsi_detail::Grid grid =
        Tesseroid::tli_mmsi_detail::make_grid(
            MODEL_LONS, MODEL_LONE, MODEL_LATS, MODEL_LATE, MODEL_DLON,
            MODEL_DLAT, model_count);

    // ---------- Observation grid (lon-major, matching TLI-MMSI storage) ----------
    int sampled_count = 0;
    for (int i = 0; i < OBS_COUNT; ++i)
        if (i % stride == 0)
            ++sampled_count;

    std::vector<double> obs_lon(static_cast<std::size_t>(sampled_count));
    std::vector<double> obs_lat(static_cast<std::size_t>(sampled_count));

    if (rank == 0) {
        int p = 0;
        for (int ilon = 0; ilon < OBS_NLON; ++ilon) {
            const double lon = OBS_LON_FIRST + ilon * 1.0;
            for (int ilat = 0; ilat < OBS_NLAT; ++ilat) {
                if (p % stride != 0) {
                    ++p;
                    continue;
                }
                const std::size_t q = static_cast<std::size_t>(p / stride);
                obs_lon[q] = lon;
                obs_lat[q] = OBS_LAT_FIRST + ilat * 1.0;
                ++p;
            }
        }
    }
    MPI_Bcast(obs_lon.data(), sampled_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(obs_lat.data(), sampled_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    const double heights[3] = {R_E + 1.0, R_E - 31.0, R_E - 15.0};
    const char* test_name[3] = {"Test5_outside_plus1km",
                                "Test6_inside_minus31km",
                                "Test7_within_minus15km"};
    const char* dat_suffix[3] = {"1km", "-31km", "-15km"};

    // Distribute sampled points across MPI processes.
    const int local_count_int = local_count_int_for_rank(
        rank, sampled_count, processes);
    const int offset = rank * (sampled_count / processes) +
                       std::min(rank, sampled_count % processes);
    const std::size_t local_count = static_cast<std::size_t>(local_count_int);

    // Per-test, per-component local statistics.
    std::array<std::vector<Stats>, 3> local_stats;
    for (auto& test : local_stats)
        test.resize(20);

    const auto start = std::chrono::steady_clock::now();
    char buffer[1024];

    // Strings accumulate on rank 0 for the summary CSV.
    std::string summary_csv;

    for (int t = first_test; t <= last_test; ++t) {
        const double r = heights[t];
        std::snprintf(buffer, sizeof(buffer),
                      "Test %d (%s): r=%.4f km, %d points, %d MPI processes",
                      t + 5, test_name[t], r, sampled_count, processes);
        if (rank == 0) {
            std::printf("[%8.1f s] %s\n",
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - start).count(),
                        buffer);
            std::fflush(stdout);
        }
        const auto test_start = std::chrono::steady_clock::now();

        // Per-column arrays in the manuscript Work layout:
        // lon[0..n), lat, r, then the 20 components in Work order.
        std::vector<std::vector<double>> columns(23);
        for (auto& column : columns)
            column.reserve(local_count);
        for (std::size_t i = 0; i < local_count; ++i) {
            const int global = offset + static_cast<int>(i);
            const double lat = obs_lat[static_cast<std::size_t>(global)];
            const double lon = obs_lon[static_cast<std::size_t>(global)];

            const Tesseroid::tli_mmsi_detail::Values values =
                Tesseroid::tli_mmsi_detail::point_estimate(
                    grid, model_r2, model_r1, model_rho, lat, lon, r, 0.0,
                    0.0, 20);

            columns[0].push_back(lon);
            columns[1].push_back(lat);
            columns[2].push_back(r);
            for (int k = 0; k < 20; ++k)
                columns[static_cast<std::size_t>(3 + k)].push_back(
                    values[static_cast<std::size_t>(work_order[k])]);

            const ShellReference ref_s = shell_reference(r);
            const std::array<double, 20> ref = analytical(ref_s);
            for (int k = 0; k < 20; ++k)
                accumulate(local_stats[t][static_cast<std::size_t>(k)],
                           values[static_cast<std::size_t>(k)],
                           ref[static_cast<std::size_t>(k)]);
        }

        // Each process writes its own column-array slice; rank 0 combines.
        {
            char slice_name[256];
            std::snprintf(slice_name, sizeof(slice_name),
                          "results/slice_%s_rank%02d.dat", test_name[t], rank);
            std::FILE* slice = std::fopen(slice_name, "wb");
            if (slice == nullptr) {
                std::printf("Warning: could not write local slice %s\n",
                            slice_name);
            } else {
                bool ok = std::fwrite(&local_count_int, sizeof(int), 1,
                                      slice) == 1;
                for (std::size_t col = 0; col < columns.size(); ++col)
                    ok = ok && std::fwrite(columns[col].data(), sizeof(double),
                                           columns[col].size(), slice) ==
                                  columns[col].size();
                std::fclose(slice);
                if (!ok) {
                    std::fprintf(stderr, "Error: failed writing %s\n",
                                 slice_name);
                    MPI_Abort(MPI_COMM_WORLD, 4);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Rank 0 merges the per-rank column slices into the Work-format
        // .dat file: int32 count, then 23 consecutive double arrays.
        if (rank == 0) {
            _mkdir("results");
            char out_path[256];
            std::snprintf(out_path, sizeof(out_path),
                          "results/ResultsPan_10x10Split_%s.dat",
                          dat_suffix[t]);
            std::FILE* out = std::fopen(out_path, "wb");
            if (out == nullptr) {
                std::printf("Warning: could not open %s\n", out_path);
            } else {
                std::fwrite(&sampled_count, sizeof(int), 1, out);
                for (std::size_t col = 0; col < 23; ++col) {
                    for (int p = 0; p < processes; ++p) {
                        char slice_name[256];
                        std::snprintf(slice_name, sizeof(slice_name),
                                      "results/slice_%s_rank%02d.dat",
                                      test_name[t], p);
                        std::FILE* slice = std::fopen(slice_name, "rb");
                        if (slice == nullptr) {
                            std::printf("Warning: missing slice %s\n",
                                        slice_name);
                            continue;
                        }
                        int slice_count = 0;
                        if (std::fread(&slice_count, sizeof(int), 1, slice) != 1 ||
                            slice_count != local_count_int_for_rank(p, sampled_count,
                                                                     processes)) {
                            std::fprintf(stderr, "Error: invalid slice %s\n",
                                         slice_name);
                            std::fclose(slice);
                            MPI_Abort(MPI_COMM_WORLD, 5);
                        }
                        // Skip the columns before the one we want:
                        // int32 header + col columns of slice_count doubles.
                        std::fseek(slice,
                                   static_cast<long>(sizeof(int) +
                                                     col * slice_count *
                                                         sizeof(double)),
                                   SEEK_SET);
                        std::vector<double> column(
                            static_cast<std::size_t>(slice_count));
                        if (std::fread(column.data(), sizeof(double),
                                       column.size(), slice) != column.size() ||
                            std::fwrite(column.data(), sizeof(double),
                                        column.size(), out) != column.size()) {
                            std::fprintf(stderr, "Error: corrupt slice/output while merging %s\n",
                                         slice_name);
                            std::fclose(slice);
                            std::fclose(out);
                            MPI_Abort(MPI_COMM_WORLD, 6);
                        }
                        std::fclose(slice);
                    }
                }
                std::fclose(out);
                std::printf("Wrote %s (%d points x 23 columns)\n", out_path,
                            sampled_count);
                std::fflush(stdout);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        {
            char slice_name[256];
            std::snprintf(slice_name, sizeof(slice_name),
                          "results/slice_%s_rank%02d.dat", test_name[t], rank);
            std::remove(slice_name);
        }

        // Reduce statistics over all processes.
        std::array<double, 20> global_min{}, global_max{}, global_sum{},
            global_sumsq{};
        std::array<unsigned long long, 20> global_count{}, global_nan{};
        for (int k = 0; k < 20; ++k) {
            MPI_Allreduce(&local_stats[t][static_cast<std::size_t>(k)].min,
                          &global_min[static_cast<std::size_t>(k)], 1,
                          MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&local_stats[t][static_cast<std::size_t>(k)].max,
                          &global_max[static_cast<std::size_t>(k)], 1,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&local_stats[t][static_cast<std::size_t>(k)].sum,
                          &global_sum[static_cast<std::size_t>(k)], 1,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_stats[t][static_cast<std::size_t>(k)].sumsq,
                          &global_sumsq[static_cast<std::size_t>(k)], 1,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            unsigned long long local_count_value =
                static_cast<unsigned long long>(
                    local_stats[t][static_cast<std::size_t>(k)].count);
            MPI_Allreduce(&local_count_value,
                          &global_count[static_cast<std::size_t>(k)], 1,
                          MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            unsigned long long local_nan_value =
                static_cast<unsigned long long>(
                    local_stats[t][static_cast<std::size_t>(k)].nan_count);
            MPI_Allreduce(&local_nan_value,
                          &global_nan[static_cast<std::size_t>(k)], 1,
                          MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        }

        if (rank == 0) {
            const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - test_start).count();
            std::printf("[%8.1f s] Test %s done: %zu points, elapsed %.1f s, "
                        "%.4f ms/point\n",
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - start).count(),
                        test_name[t], global_count[0], seconds,
                        1000.0 * seconds /
                            std::max<unsigned long long>(1, global_count[0]));
            std::printf("%-7s %-5s %8s %13s %13s %13s %13s\n", "comp",
                        "order", "count", "min_abs", "max_abs", "mean_abs",
                        "rmse_abs");
            for (int k = 0; k < 20; ++k) {
                const unsigned long long n =
                    global_count[static_cast<std::size_t>(k)];
                const double mean =
                    n > 0
                        ? global_sum[static_cast<std::size_t>(k)] /
                              static_cast<double>(n)
                        : 0.0;
                const double rmse =
                    n > 0
                        ? std::sqrt(global_sumsq[static_cast<std::size_t>(k)] /
                                    static_cast<double>(n))
                        : 0.0;
                std::printf("%-7s %-5d %8llu %13.6e %13.6e %13.6e %13.6e\n",
                            component_name(k), component_order(k), n,
                            global_min[static_cast<std::size_t>(k)],
                            global_max[static_cast<std::size_t>(k)], mean,
                            rmse);
            }
            std::printf("\n");
            std::fflush(stdout);
        }

        // Collect this test's summary rows on rank 0 for the CSV file.
        if (rank == 0) {
            for (int k = 0; k < 20; ++k) {
                const unsigned long long n =
                    global_count[static_cast<std::size_t>(k)];
                const double mean =
                    n > 0
                        ? global_sum[static_cast<std::size_t>(k)] /
                              static_cast<double>(n)
                        : 0.0;
                const double rmse =
                    n > 0
                        ? std::sqrt(global_sumsq[static_cast<std::size_t>(k)] /
                                    static_cast<double>(n))
                        : 0.0;
                std::snprintf(buffer, sizeof(buffer),
                              "%s,%.4f,%s,%d,%llu,%llu,"
                              "%.6e,%.6e,%.6e,%.6e\n",
                              test_name[t], r, component_name(k),
                              component_order(k), n,
                              global_nan[static_cast<std::size_t>(k)],
                              global_min[static_cast<std::size_t>(k)],
                              global_max[static_cast<std::size_t>(k)], mean,
                              rmse);
                summary_csv += buffer;
            }
        }
        if (t == 2 && rank == 0) {
            _mkdir("results");
            std::FILE* file = std::fopen("results/shell_summary.csv", "w");
            if (file != nullptr) {
                std::fputs("test,r_km,component,order,count,nan,"
                           "min_abs,max_abs,mean_abs,rmse_abs\n",
                           file);
                std::fputs(summary_csv.c_str(), file);
                std::fclose(file);
                std::printf("Wrote results/shell_summary.csv\n");
            }
        }
    }

    if (rank == 0) {
        const double total = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        std::printf("Total elapsed: %.1f s\n", total);
    }

    MPI_Finalize();
    return 0;
}
