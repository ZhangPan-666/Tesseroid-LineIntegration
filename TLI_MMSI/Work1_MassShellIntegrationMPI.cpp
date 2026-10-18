// Run Work1 mass-shell integration with MPI.

#include "_tesseroid_MSI.h"
#include <algorithm>
#include <chrono>

class Work1MassTimer {
  public:
    void start() { start_ = Clock::now(); }
    double seconds() const {
        return std::chrono::duration<double>(Clock::now() - start_).count();
    }

  private:
    using Clock = std::chrono::high_resolution_clock;
    std::chrono::time_point<Clock> start_;
};

int main(int argc, char** argv) {
    int rank, nprocs;
    int NumLon, NumLat, NumCoordinates;
    long long NumModel;

    double AbsTol, RelTol;
    double Lons, Lone, Lats, Late, dLon, dLat;

    double *R1, *R2, *rho;
    double *Longitude, *Latitude, *Radius;
    double *V, *Vx, *Vy, *Vz, *Vxx, *Vxy, *Vyy, *Vzx, *Vzy, *Vzz, *Vxxx, *Vxxy, *Vxxz, *Vxyz, *Vyyx, *Vyyy, *Vyyz,
        *Vzzx, *Vzzy, *Vzzz;

    char ParametersFile[1024] = "MSIG.ForPar";
    char ResultsFile[1024] = "ResultG.dat";

    FILE* fp_Parameters;
    FILE* fp_Results;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    R1 = NULL;
    R2 = NULL;
    rho = NULL;
    Longitude = NULL;
    Latitude = NULL;
    Radius = NULL;
    V = NULL;
    Vx = Vy = Vz = NULL;
    Vxx = Vxy = Vyy = Vzx = Vzy = Vzz = NULL;
    Vxxx = Vxxy = Vxxz = Vxyz = Vyyx = Vyyy = Vyyz = Vzzx = Vzzy = Vzzz = NULL;

    int ReadOK = 1;

    if (rank == 0) {
        fp_Parameters = fopen(ParametersFile, "rb");

        if (fp_Parameters == NULL) {
            printf("\nError: failed to open parameter file.\n");
            printf("Calculation aborted.\n");
            ReadOK = 0;
        }

        if (ReadOK) {
            if (fread(&AbsTol, 8, 1, fp_Parameters) != 1 || fread(&RelTol, 8, 1, fp_Parameters) != 1 ||
                fread(&Lons, 8, 1, fp_Parameters) != 1 || fread(&Lone, 8, 1, fp_Parameters) != 1 ||
                fread(&Lats, 8, 1, fp_Parameters) != 1 || fread(&Late, 8, 1, fp_Parameters) != 1 ||
                fread(&dLon, 8, 1, fp_Parameters) != 1 || fread(&dLat, 8, 1, fp_Parameters) != 1 ||
                fread(&NumLon, 4, 1, fp_Parameters) != 1 || fread(&NumLat, 4, 1, fp_Parameters) != 1) {
                printf("\nError: failed to read parameter file header.\n");
                ReadOK = 0;
            }
        }

        if (ReadOK) {
            NumModel = static_cast<long long>(NumLon) * static_cast<long long>(NumLat);

            R1 = (double*)malloc(sizeof(double) * NumModel);
            R2 = (double*)malloc(sizeof(double) * NumModel);
            rho = (double*)malloc(sizeof(double) * NumModel);

            if (R1 == NULL || R2 == NULL || rho == NULL) {
                printf("\nError: failed to allocate model memory.\n");
                ReadOK = 0;
            }
        }

        if (ReadOK) {
            if (fread(R2, 8, NumModel, fp_Parameters) != (size_t)NumModel ||
                fread(R1, 8, NumModel, fp_Parameters) != (size_t)NumModel ||
                fread(rho, 8, NumModel, fp_Parameters) != (size_t)NumModel ||
                fread(&NumCoordinates, 4, 1, fp_Parameters) != 1) {
                printf("\nError: failed to read model data.\n");
                ReadOK = 0;
            }
        }

        if (ReadOK) {
            Longitude = (double*)malloc(sizeof(double) * NumCoordinates);
            Latitude = (double*)malloc(sizeof(double) * NumCoordinates);
            Radius = (double*)malloc(sizeof(double) * NumCoordinates);

            if (Longitude == NULL || Latitude == NULL || Radius == NULL) {
                printf("\nError: failed to allocate coordinate memory.\n");
                ReadOK = 0;
            }
        }

        if (ReadOK) {
            if (fread(Longitude, 8, NumCoordinates, fp_Parameters) != (size_t)NumCoordinates ||
                fread(Latitude, 8, NumCoordinates, fp_Parameters) != (size_t)NumCoordinates ||
                fread(Radius, 8, NumCoordinates, fp_Parameters) != (size_t)NumCoordinates) {
                printf("\nError: failed to read coordinate data.\n");
                ReadOK = 0;
            }
        }

        if (fp_Parameters != NULL) {
            fclose(fp_Parameters);
        }
    }

    MPI_Bcast(&ReadOK, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (!ReadOK) {
        free(R1);
        free(R2);
        free(rho);
        free(Longitude);
        free(Latitude);
        free(Radius);
        MPI_Finalize();
        return 0;
    }

    MPI_Bcast(&AbsTol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&RelTol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&Lons, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Lone, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Lats, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Late, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dLon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dLat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&NumLon, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NumLat, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NumCoordinates, 1, MPI_INT, 0, MPI_COMM_WORLD);

    NumModel = static_cast<long long>(NumLon) * static_cast<long long>(NumLat);

    if (rank != 0) {
        R1 = (double*)malloc(sizeof(double) * NumModel);
        R2 = (double*)malloc(sizeof(double) * NumModel);
        rho = (double*)malloc(sizeof(double) * NumModel);

        Longitude = (double*)malloc(sizeof(double) * NumCoordinates);
        Latitude = (double*)malloc(sizeof(double) * NumCoordinates);
        Radius = (double*)malloc(sizeof(double) * NumCoordinates);
    }

    if (R1 == NULL || R2 == NULL || rho == NULL || Longitude == NULL || Latitude == NULL || Radius == NULL) {
        if (rank == 0) {
            printf("\nError: failed to allocate memory after MPI broadcast.\n");
        }
        free(R1);
        free(R2);
        free(rho);
        free(Longitude);
        free(Latitude);
        free(Radius);
        MPI_Finalize();
        return 0;
    }

    V = (double*)malloc(sizeof(double) * NumCoordinates);

    Vx = (double*)malloc(sizeof(double) * NumCoordinates);
    Vy = (double*)malloc(sizeof(double) * NumCoordinates);
    Vz = (double*)malloc(sizeof(double) * NumCoordinates);

    Vxx = (double*)malloc(sizeof(double) * NumCoordinates);
    Vxy = (double*)malloc(sizeof(double) * NumCoordinates);
    Vyy = (double*)malloc(sizeof(double) * NumCoordinates);
    Vzx = (double*)malloc(sizeof(double) * NumCoordinates);
    Vzy = (double*)malloc(sizeof(double) * NumCoordinates);
    Vzz = (double*)malloc(sizeof(double) * NumCoordinates);

    Vxxx = (double*)malloc(sizeof(double) * NumCoordinates);
    Vxxy = (double*)malloc(sizeof(double) * NumCoordinates);
    Vxxz = (double*)malloc(sizeof(double) * NumCoordinates);
    Vxyz = (double*)malloc(sizeof(double) * NumCoordinates);
    Vyyx = (double*)malloc(sizeof(double) * NumCoordinates);
    Vyyy = (double*)malloc(sizeof(double) * NumCoordinates);
    Vyyz = (double*)malloc(sizeof(double) * NumCoordinates);
    Vzzx = (double*)malloc(sizeof(double) * NumCoordinates);
    Vzzy = (double*)malloc(sizeof(double) * NumCoordinates);
    Vzzz = (double*)malloc(sizeof(double) * NumCoordinates);

    if (V == NULL || Vx == NULL || Vy == NULL || Vz == NULL || Vxx == NULL || Vxy == NULL || Vyy == NULL ||
        Vzx == NULL || Vzy == NULL || Vzz == NULL || Vxxx == NULL || Vxxy == NULL || Vxxz == NULL || Vxyz == NULL ||
        Vyyx == NULL || Vyyy == NULL || Vyyz == NULL || Vzzx == NULL || Vzzy == NULL || Vzzz == NULL) {
        if (rank == 0) {
            printf("\nError: failed to allocate result memory.\n");
        }
        free(R1);
        free(R2);
        free(rho);
        free(Longitude);
        free(Latitude);
        free(Radius);
        free(V);
        free(Vx);
        free(Vy);
        free(Vz);
        free(Vxx);
        free(Vxy);
        free(Vyy);
        free(Vzx);
        free(Vzy);
        free(Vzz);
        free(Vxxx);
        free(Vxxy);
        free(Vxxz);
        free(Vxyz);
        free(Vyyx);
        free(Vyyy);
        free(Vyyz);
        free(Vzzx);
        free(Vzzy);
        free(Vzzz);
        MPI_Finalize();
        return 0;
    }

    if (rank == 0) {
        const auto r2Range = std::minmax_element(R2, R2 + NumModel);
        const auto r1Range = std::minmax_element(R1, R1 + NumModel);
        const auto rhoRange = std::minmax_element(rho, rho + NumModel);
        const auto lonRange = std::minmax_element(Longitude, Longitude + NumCoordinates);
        const auto latRange = std::minmax_element(Latitude, Latitude + NumCoordinates);
        const auto radRange = std::minmax_element(Radius, Radius + NumCoordinates);
        printf("\n========== Work1 mass-shell parameters ==========\n");
        printf("Files: %s -> %s | MPI processes: %d\n", ParametersFile, ResultsFile, nprocs);
        printf("Grid: lon [%.15g, %.15g], step %.15g deg; lat [%.15g, %.15g], step %.15g deg\n", Lons, Lone, dLon, Lats,
               Late, dLat);
        printf("Grid size: %d x %d = %lld nodes\n", NumLon, NumLat, NumModel);
        printf("R2 range: [%.15g, %.15g]; R1 range: [%.15g, %.15g]\n", *r2Range.first, *r2Range.second, *r1Range.first,
               *r1Range.second);
        printf("Density range: [%.15g, %.15g]\n", *rhoRange.first, *rhoRange.second);
        printf("Observations: %d; lon [%.15g, %.15g] deg; lat [%.15g, %.15g] deg; radius [%.15g, %.15g]\n",
               NumCoordinates, *lonRange.first, *lonRange.second, *latRange.first, *latRange.second, *radRange.first,
               *radRange.second);
        printf("Tolerance: abs %.6e, rel %.6e\n", AbsTol, RelTol);
        printf("Components: V, first-, second-, and third-order derivatives\n");
        printf("==================================================\n\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    Work1MassTimer timer;
    timer.start();
    Tesseroid::MassShellIntegrationMPI(Lons, Lone, Lats, Late, dLon, dLat, R2, R1, rho, NumCoordinates, Latitude,
                                       Longitude, Radius, AbsTol, RelTol, V, Vx, Vy, Vz, Vxx, Vxy, Vyy, Vzx, Vzy, Vzz,
                                       Vxxx, Vxxy, Vxxz, Vxyz, Vyyx, Vyyy, Vyyz, Vzzx, Vzzy, Vzzz);

    MPI_Barrier(MPI_COMM_WORLD);
    const double localElapsed = timer.seconds();
    double maximumElapsed = 0.0;
    MPI_Reduce(&localElapsed, &maximumElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
        printf("Work1 TLI-MMSI gravity elapsed time: %.6f seconds\n", maximumElapsed);

    if (rank == 0) {
        fp_Results = fopen(ResultsFile, "wb");

        if (fp_Results == NULL) {
            printf("\nError: failed to open result file.\n");
        } else {
            if (fwrite(&NumCoordinates, 4, 1, fp_Results) != 1 ||
                fwrite(Longitude, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Latitude, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Radius, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(V, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vx, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vy, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vz, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vxx, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vxy, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vyy, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vzx, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vzy, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vzz, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vxxx, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vxxy, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vxxz, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vxyz, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vyyx, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vyyy, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vyyz, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vzzx, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vzzy, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates ||
                fwrite(Vzzz, 8, NumCoordinates, fp_Results) != (size_t)NumCoordinates) {
                printf("\nError: failed to write result file.\n");
            }

            fclose(fp_Results);
        }
    }

    free(R1);
    free(R2);
    free(rho);

    free(Longitude);
    free(Latitude);
    free(Radius);

    free(V);
    free(Vx);
    free(Vy);
    free(Vz);
    free(Vxx);
    free(Vxy);
    free(Vyy);
    free(Vzx);
    free(Vzy);
    free(Vzz);
    free(Vxxx);
    free(Vxxy);
    free(Vxxz);
    free(Vxyz);
    free(Vyyx);
    free(Vyyy);
    free(Vyyz);
    free(Vzzx);
    free(Vzzy);
    free(Vzzz);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
