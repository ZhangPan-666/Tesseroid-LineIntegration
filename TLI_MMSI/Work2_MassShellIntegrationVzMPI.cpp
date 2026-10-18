// Run Work2 mass-shell integration for Vz with MPI.
// The result file keeps the Work1 mass-result binary layout.

#include "_tesseroid_MSI.h"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <limits>
#include <vector>

namespace Work2MassVz {

inline void Work2MassVzComputeMPI(double lons, double lone, double lats, double late, double dLon, double dLat,
                                  std::vector<double>& outerRadius, std::vector<double>& innerRadius,
                                  std::vector<double>& density, std::vector<double>& latitude,
                                  std::vector<double>& longitude, std::vector<double>& radius, double absoluteTolerance,
                                  double relativeTolerance, std::vector<double>& vz) {
    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const int pointCount = static_cast<int>(latitude.size());
    const int modelCount = static_cast<int>(outerRadius.size());
    MPI_Barrier(MPI_COMM_WORLD);
    const auto start = std::chrono::steady_clock::now();
    MPI_Bcast(outerRadius.data(), modelCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(innerRadius.data(), modelCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(density.data(), modelCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(latitude.data(), pointCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(longitude.data(), pointCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(radius.data(), pointCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::vector<int> counts(nprocs), displacements(nprocs);
    for (int p = 0; p < nprocs; ++p) {
        counts[p] = pointCount / nprocs + (p < pointCount % nprocs ? 1 : 0);
        displacements[p] = p == 0 ? 0 : displacements[p - 1] + counts[p - 1];
    }
    const int localCount = counts[rank];
    const int offset = displacements[rank];
    std::vector<double> localVz(localCount);
    const auto grid = Tesseroid::tli_mmsi_detail::make_grid(
        lons, lone, lats, late, dLon, dLat, outerRadius.size());

    for (int i = 0; i < localCount; ++i) {
        const int globalIndex = offset + i;
        const auto values =
            Tesseroid::tli_mmsi_detail::point_estimate(
            grid, outerRadius, innerRadius, density, latitude[globalIndex],
            longitude[globalIndex], radius[globalIndex], absoluteTolerance,
            relativeTolerance, 3, 1u << 3);
        localVz[i] = values[3];
    }

    if (rank == 0)
        vz.resize(pointCount);
    MPI_Gatherv(localVz.data(), localCount, MPI_DOUBLE,
                rank == 0 ? vz.data() : nullptr, counts.data(),
                displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    const double localSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    double maximumSeconds = 0.0;
    MPI_Reduce(&localSeconds, &maximumSeconds, 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Work2 TLI-MMSI gravity Vz elapsed time: "
                  << maximumSeconds << " seconds\n";
    }
}

}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const char* parameterFile = "MSIG.ForPar";
    const char* resultFile = "ResultG.dat";

    double absoluteTolerance = 0.0;
    double relativeTolerance = 0.0;
    double lons = 0.0;
    double lone = 0.0;
    double lats = 0.0;
    double late = 0.0;
    double dLon = 0.0;
    double dLat = 0.0;
    int numLon = 0;
    int numLat = 0;
    int pointCount = 0;
    int readOK = 1;

    std::vector<double> outerRadius;
    std::vector<double> innerRadius;
    std::vector<double> density;
    std::vector<double> longitude;
    std::vector<double> latitude;
    std::vector<double> radius;

    if (rank == 0) {
        FILE* file = std::fopen(parameterFile, "rb");
        if (file == nullptr) {
            std::cerr << "Error: failed to open " << parameterFile << "\n";
            readOK = 0;
        }

        if (readOK) {
            readOK =
                std::fread(&absoluteTolerance, sizeof(double), 1, file) == 1 &&
                std::fread(&relativeTolerance, sizeof(double), 1, file) == 1 &&
                std::fread(&lons, sizeof(double), 1, file) == 1 && std::fread(&lone, sizeof(double), 1, file) == 1 &&
                std::fread(&lats, sizeof(double), 1, file) == 1 && std::fread(&late, sizeof(double), 1, file) == 1 &&
                std::fread(&dLon, sizeof(double), 1, file) == 1 && std::fread(&dLat, sizeof(double), 1, file) == 1 &&
                std::fread(&numLon, sizeof(int), 1, file) == 1 && std::fread(&numLat, sizeof(int), 1, file) == 1;
        }

        if (readOK && (numLon < 2 || numLat < 2 || dLon <= 0.0 || dLat <= 0.0)) {
            std::cerr << "Error: invalid model dimensions or grid spacing.\n";
            readOK = 0;
        }

        if (readOK) {
            const std::size_t modelCount = static_cast<std::size_t>(numLon) * static_cast<std::size_t>(numLat);
            outerRadius.resize(modelCount);
            innerRadius.resize(modelCount);
            density.resize(modelCount);
            readOK = std::fread(outerRadius.data(), sizeof(double), modelCount, file) == modelCount &&
                     std::fread(innerRadius.data(), sizeof(double), modelCount, file) == modelCount &&
                     std::fread(density.data(), sizeof(double), modelCount, file) == modelCount &&
                     std::fread(&pointCount, sizeof(int), 1, file) == 1;
        }

        if (readOK && pointCount > 0) {
            longitude.resize(pointCount);
            latitude.resize(pointCount);
            radius.resize(pointCount);
            readOK =
                std::fread(longitude.data(), sizeof(double), pointCount, file) ==
                    static_cast<std::size_t>(pointCount) &&
                std::fread(latitude.data(), sizeof(double), pointCount, file) == static_cast<std::size_t>(pointCount) &&
                std::fread(radius.data(), sizeof(double), pointCount, file) == static_cast<std::size_t>(pointCount);
        } else if (readOK) {
            std::cerr << "Error: no observation points were supplied.\n";
            readOK = 0;
        }

        if (file != nullptr)
            std::fclose(file);
    }

    if (rank == 0 && readOK) {
        const auto outer = std::minmax_element(outerRadius.begin(), outerRadius.end());
        const auto inner = std::minmax_element(innerRadius.begin(), innerRadius.end());
        const auto rho = std::minmax_element(density.begin(), density.end());
        const auto obsLon = std::minmax_element(longitude.begin(), longitude.end());
        const auto obsLat = std::minmax_element(latitude.begin(), latitude.end());
        const auto obsRad = std::minmax_element(radius.begin(), radius.end());
        std::cout << "\n========== Work2 mass-shell parameters ==========\n"
                  << "Files: " << parameterFile << " -> " << resultFile << " | MPI processes: " << nprocs << '\n'
                  << "Grid: lon [" << lons << ", " << lone << "], step " << dLon << " deg; lat [" << lats << ", "
                  << late << "], step " << dLat << " deg\n"
                  << "Grid size: " << numLon << " x " << numLat << " = " << outerRadius.size() << " nodes\n"
                  << "Outer radius: [" << *outer.first << ", " << *outer.second << "]; inner radius: [" << *inner.first
                  << ", " << *inner.second << "]\n"
                  << "Density range: [" << *rho.first << ", " << *rho.second << "]\n"
                  << "Observations: " << pointCount << "; lon [" << *obsLon.first << ", " << *obsLon.second
                  << "] deg; lat [" << *obsLat.first << ", " << *obsLat.second << "] deg; radius [" << *obsRad.first
                  << ", " << *obsRad.second << "]\n"
                  << "Tolerance: abs " << absoluteTolerance << ", rel " << relativeTolerance << '\n'
                  << "Components: Vz only; other Work1-layout fields are NaN\n"
                  << "=================================================\n\n"
                  << std::flush;
        std::cout << "Starting Work2 mass-shell Vz integration...\n\n";
    }

    MPI_Bcast(&readOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!readOK) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    MPI_Bcast(&absoluteTolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&relativeTolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lons, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lone, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lats, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&late, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dLon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dLat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numLon, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numLat, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&pointCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        const std::size_t modelCount = static_cast<std::size_t>(numLon) * static_cast<std::size_t>(numLat);
        outerRadius.resize(modelCount);
        innerRadius.resize(modelCount);
        density.resize(modelCount);
        longitude.resize(pointCount);
        latitude.resize(pointCount);
        radius.resize(pointCount);
    }

    std::vector<double> vz;
    Work2MassVz::Work2MassVzComputeMPI(lons, lone, lats, late, dLon, dLat, outerRadius, innerRadius, density, latitude,
                                       longitude, radius, absoluteTolerance, relativeTolerance, vz);

    int writeOK = 1;
    if (rank == 0) {
        const std::vector<double> notComputed(pointCount, std::numeric_limits<double>::quiet_NaN());
        auto writeArray = [](FILE* file, const std::vector<double>& values, int count) {
            return std::fwrite(values.data(), sizeof(double), count, file) == static_cast<std::size_t>(count);
        };

        FILE* file = std::fopen(resultFile, "wb");
        if (file == nullptr) {
            std::cerr << "Error: failed to open " << resultFile << "\n";
            writeOK = 0;
        }
        if (writeOK) {
            writeOK = std::fwrite(&pointCount, sizeof(int), 1, file) == 1 && writeArray(file, longitude, pointCount) &&
                      writeArray(file, latitude, pointCount) && writeArray(file, radius, pointCount) &&

                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, vz, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount) &&
                      writeArray(file, notComputed, pointCount) && writeArray(file, notComputed, pointCount);
        }
        if (file != nullptr)
            std::fclose(file);
    }

    MPI_Bcast(&writeOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return writeOK ? EXIT_SUCCESS : EXIT_FAILURE;
}
