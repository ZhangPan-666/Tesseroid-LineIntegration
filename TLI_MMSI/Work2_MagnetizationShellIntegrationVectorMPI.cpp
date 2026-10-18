// Run Work2 magnetization-shell integration for Vx, Vy, and Vz with MPI.
// The result file keeps the Work1 magnetic-result binary layout.

#include "_tesseroid_MSI.h"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <vector>

namespace Work2MagVector {

inline void Work2MagVectorComputeMPI(double lons, double lone, double lats, double late, double dLon, double dLat,
                                     std::vector<double>& outerRadius, std::vector<double>& innerRadius,
                                     std::vector<double>& magnetizationX, std::vector<double>& magnetizationY,
                                     std::vector<double>& magnetizationZ, std::vector<double>& latitude,
                                     std::vector<double>& longitude, std::vector<double>& radius,
                                     double absoluteTolerance, double relativeTolerance, std::vector<double>& vx,
                                     std::vector<double>& vy, std::vector<double>& vz) {
    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const int pointCount = static_cast<int>(latitude.size());
    const int modelCount = static_cast<int>(outerRadius.size());
    MPI_Barrier(MPI_COMM_WORLD);
    const auto start = std::chrono::steady_clock::now();
    for (double* source : {outerRadius.data(), innerRadius.data(),
                           magnetizationX.data(), magnetizationY.data(),
                           magnetizationZ.data()})
        MPI_Bcast(source, modelCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
    std::vector<double> localVx(localCount), localVy(localCount),
        localVz(localCount);
    const auto grid = Tesseroid::tli_mmsi_detail::make_grid(
        lons, lone, lats, late, dLon, dLat, outerRadius.size());

    // Reuse the transformed magnetization arrays for every observation handled
    // by this MPI rank.  The former implementation allocated all three arrays
    // inside Tesseroid_MagneticPointEstimate for every observation.
    std::vector<double> mxLocal(modelCount), myLocal(modelCount),
        mzLocal(modelCount);
    std::vector<double> sourceLonCos(grid.nlon), sourceLonSin(grid.nlon);
    std::vector<double> sourceLatCos(grid.nlat), sourceLatSin(grid.nlat);
    for (std::size_t i = 0; i < grid.nlon; ++i) {
        const double sourceLongitude =
            (lons + static_cast<double>(i) * dLon) *
            Tesseroid::tli_mmsi_detail::degree;
        sourceLonCos[i] = std::cos(sourceLongitude);
        sourceLonSin[i] = std::sin(sourceLongitude);
    }
    for (std::size_t j = 0; j < grid.nlat; ++j) {
        const double sourceLatitude =
            (lats + static_cast<double>(j) * dLat) *
            Tesseroid::tli_mmsi_detail::degree;
        sourceLatCos[j] = std::cos(sourceLatitude);
        sourceLatSin[j] = std::sin(sourceLatitude);
    }

    for (int i = 0; i < localCount; ++i) {
        const int globalIndex = offset + i;
        const double observationLatitude =
            latitude[globalIndex] * Tesseroid::tli_mmsi_detail::degree;
        const double observationLongitude =
            longitude[globalIndex] * Tesseroid::tli_mmsi_detail::degree;
        const double co = std::cos(observationLatitude);
        const double so = std::sin(observationLatitude);
        const double observationLonCos = std::cos(observationLongitude);
        const double observationLonSin = std::sin(observationLongitude);

        for (std::size_t ilon = 0; ilon < grid.nlon; ++ilon) {
            const double cb = sourceLonCos[ilon] * observationLonCos +
                              sourceLonSin[ilon] * observationLonSin;
            const double sb = sourceLonSin[ilon] * observationLonCos -
                              sourceLonCos[ilon] * observationLonSin;
            for (std::size_t ilat = 0; ilat < grid.nlat; ++ilat) {
                const std::size_t p = ilon * grid.nlat + ilat;
                const double cs = sourceLatCos[ilat];
                const double ss = sourceLatSin[ilat];
                mxLocal[p] = (cb * so * ss + co * cs) * magnetizationX[p] +
                             sb * so * magnetizationY[p] -
                             (cb * so * cs - co * ss) * magnetizationZ[p];
                myLocal[p] = -sb * ss * magnetizationX[p] +
                             cb * magnetizationY[p] +
                             sb * cs * magnetizationZ[p];
                mzLocal[p] = -(cb * co * ss - so * cs) * magnetizationX[p] -
                             sb * co * magnetizationY[p] +
                             (cb * co * cs + so * ss) * magnetizationZ[p];
            }
        }

        auto estimateComponent = [&](const std::vector<double>& source,
                                     std::size_t component) {
            const auto values = Tesseroid::tli_mmsi_detail::point_estimate(
                grid, outerRadius, innerRadius, source,
                latitude[globalIndex], longitude[globalIndex],
                radius[globalIndex], absoluteTolerance, relativeTolerance,
                component, 1u << component);
            return values[component];
        };
        const double gx4 = estimateComponent(mxLocal, 4);
        const double gx5 = estimateComponent(mxLocal, 5);
        const double gx6 = estimateComponent(mxLocal, 6);
        const double gy5 = estimateComponent(myLocal, 5);
        const double gy7 = estimateComponent(myLocal, 7);
        const double gy8 = estimateComponent(myLocal, 8);
        const double gz6 = estimateComponent(mzLocal, 6);
        const double gz8 = estimateComponent(mzLocal, 8);
        const double gz9 = estimateComponent(mzLocal, 9);
        localVx[i] = gx4 + gy5 + gz6;
        localVy[i] = gx5 + gy7 + gz8;
        localVz[i] = gx6 + gy8 + gz9;
    }

    if (rank == 0) {
        vx.resize(pointCount);
        vy.resize(pointCount);
        vz.resize(pointCount);
    }
    MPI_Gatherv(localVx.data(), localCount, MPI_DOUBLE,
                rank == 0 ? vx.data() : nullptr, counts.data(),
                displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(localVy.data(), localCount, MPI_DOUBLE,
                rank == 0 ? vy.data() : nullptr, counts.data(),
                displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
        std::cout << "Work2 TLI-MMSI magnetic vector elapsed time: "
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

    const char* parameterFile = "MSIM.ForPar";
    const char* resultFile = "ResultM.dat";

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
    std::vector<double> magnetizationX;
    std::vector<double> magnetizationY;
    std::vector<double> magnetizationZ;
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
            magnetizationX.resize(modelCount);
            magnetizationY.resize(modelCount);
            magnetizationZ.resize(modelCount);
            readOK = std::fread(outerRadius.data(), sizeof(double), modelCount, file) == modelCount &&
                     std::fread(innerRadius.data(), sizeof(double), modelCount, file) == modelCount &&
                     std::fread(magnetizationX.data(), sizeof(double), modelCount, file) == modelCount &&
                     std::fread(magnetizationY.data(), sizeof(double), modelCount, file) == modelCount &&
                     std::fread(magnetizationZ.data(), sizeof(double), modelCount, file) == modelCount &&
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
        const auto mx = std::minmax_element(magnetizationX.begin(), magnetizationX.end());
        const auto my = std::minmax_element(magnetizationY.begin(), magnetizationY.end());
        const auto mz = std::minmax_element(magnetizationZ.begin(), magnetizationZ.end());
        const auto obsLon = std::minmax_element(longitude.begin(), longitude.end());
        const auto obsLat = std::minmax_element(latitude.begin(), latitude.end());
        const auto obsRad = std::minmax_element(radius.begin(), radius.end());
        std::cout << "\n======= Work2 magnetization-shell parameters =======\n"
                  << "Files: " << parameterFile << " -> " << resultFile << " | MPI processes: " << nprocs << '\n'
                  << "Grid: lon [" << lons << ", " << lone << "], step " << dLon << " deg; lat [" << lats << ", "
                  << late << "], step " << dLat << " deg\n"
                  << "Grid size: " << numLon << " x " << numLat << " = " << outerRadius.size() << " nodes\n"
                  << "Outer radius: [" << *outer.first << ", " << *outer.second << "]; inner radius: [" << *inner.first
                  << ", " << *inner.second << "]\n"
                  << "Magnetization: Mx [" << *mx.first << ", " << *mx.second << "], My [" << *my.first << ", "
                  << *my.second << "], Mz [" << *mz.first << ", " << *mz.second << "]\n"
                  << "Observations: " << pointCount << "; lon [" << *obsLon.first << ", " << *obsLon.second
                  << "] deg; lat [" << *obsLat.first << ", " << *obsLat.second << "] deg; radius [" << *obsRad.first
                  << ", " << *obsRad.second << "]\n"
                  << "Tolerance: abs " << absoluteTolerance << ", rel " << relativeTolerance << '\n'
                  << "Components: Vx, Vy, Vz only; other Work1-layout fields are NaN\n"
                  << "====================================================\n\n"
                  << std::flush;
        std::cout << "Starting Work2 magnetization-shell vector integration...\n\n";
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
        magnetizationX.resize(modelCount);
        magnetizationY.resize(modelCount);
        magnetizationZ.resize(modelCount);
        longitude.resize(pointCount);
        latitude.resize(pointCount);
        radius.resize(pointCount);
    }

    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    Work2MagVector::Work2MagVectorComputeMPI(lons, lone, lats, late, dLon, dLat, outerRadius, innerRadius,
                                             magnetizationX, magnetizationY, magnetizationZ, latitude, longitude,
                                             radius, absoluteTolerance, relativeTolerance, vx, vy, vz);

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

                      writeArray(file, notComputed, pointCount) && writeArray(file, vx, pointCount) &&
                      writeArray(file, vy, pointCount) && writeArray(file, vz, pointCount) &&
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
