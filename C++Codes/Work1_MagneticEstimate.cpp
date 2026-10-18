// Run Work1 tesseroid magnetic forward modelling in serial mode.

#include <cstdio>
#include <cstdlib>

#include "_TLI_estimate.h"

int main(int argc, char** argv) {
    int i;
    int NumPermeabilityModel, NumCoordinates;

    double AbsTol, RelTol;

    double *Fai1, *Fai2, *Lamda1, *Lamda2, *R1, *R2, *Mx, *My, *Mz;
    double *Longitude, *Latitude, *Radius;
    double *V, *Vx, *Vy, *Vz, *Vxx, *Vxy, *Vyy, *Vzx, *Vzy, *Vzz;

    char ParametersFile[1024] = "TFM.ForPar";
    char ResultsFile[1024] = "Result_TFM.dat";

    FILE* fp_Parameters;
    FILE* fp_Results;

    fp_Parameters = fopen(ParametersFile, "rb");

    if (fp_Parameters == NULL) {
        printf("\nError: failed to open parameter file.\n");
        printf("Calculation aborted.\n");
        return 0;
    }

    fread(&AbsTol, 8, 1, fp_Parameters);
    fread(&RelTol, 8, 1, fp_Parameters);

    fread(&NumPermeabilityModel, 4, 1, fp_Parameters);

    Fai1 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
    Fai2 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
    Lamda1 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
    Lamda2 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
    R1 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
    R2 = (double*)malloc(sizeof(double) * NumPermeabilityModel);
    Mx = (double*)malloc(sizeof(double) * NumPermeabilityModel);
    My = (double*)malloc(sizeof(double) * NumPermeabilityModel);
    Mz = (double*)malloc(sizeof(double) * NumPermeabilityModel);

    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&Fai1[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&Fai2[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&Lamda1[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&Lamda2[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&R1[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&R2[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&Mx[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&My[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumPermeabilityModel; i++) {
        fread(&Mz[i], 8, 1, fp_Parameters);
    }

    fread(&NumCoordinates, 4, 1, fp_Parameters);

    Longitude = (double*)malloc(sizeof(double) * NumCoordinates);
    Latitude = (double*)malloc(sizeof(double) * NumCoordinates);
    Radius = (double*)malloc(sizeof(double) * NumCoordinates);

    for (i = 0; i < NumCoordinates; i++) {
        fread(&Longitude[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumCoordinates; i++) {
        fread(&Latitude[i], 8, 1, fp_Parameters);
    }
    for (i = 0; i < NumCoordinates; i++) {
        fread(&Radius[i], 8, 1, fp_Parameters);
    }

    fclose(fp_Parameters);

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

    for (i = 0; i < NumCoordinates; i++) {

        V[i] = 0.0;

        Vx[i] = 0.0;
        Vy[i] = 0.0;
        Vz[i] = 0.0;

        Vxx[i] = 0.0;
        Vxy[i] = 0.0;
        Vyy[i] = 0.0;
        Vzx[i] = 0.0;
        Vzy[i] = 0.0;
        Vzz[i] = 0.0;
    }

    Tesseroid::Tesseroid_MagneticEstimate(NumPermeabilityModel, Fai1, Fai2, Lamda1, Lamda2, R1, R2, Mx, My, Mz,
                                          NumCoordinates, Latitude, Longitude, Radius, AbsTol, RelTol, V, Vx, Vy, Vz,
                                          Vxx, Vxy, Vyy, Vzx, Vzy, Vzz);

    fp_Results = fopen(ResultsFile, "wb");

    fwrite(&NumCoordinates, 4, 1, fp_Results);

    fwrite(Longitude, 8, NumCoordinates, fp_Results);
    fwrite(Latitude, 8, NumCoordinates, fp_Results);
    fwrite(Radius, 8, NumCoordinates, fp_Results);

    fwrite(V, 8, NumCoordinates, fp_Results);
    fwrite(Vx, 8, NumCoordinates, fp_Results);
    fwrite(Vy, 8, NumCoordinates, fp_Results);
    fwrite(Vz, 8, NumCoordinates, fp_Results);
    fwrite(Vxx, 8, NumCoordinates, fp_Results);
    fwrite(Vxy, 8, NumCoordinates, fp_Results);
    fwrite(Vyy, 8, NumCoordinates, fp_Results);
    fwrite(Vzx, 8, NumCoordinates, fp_Results);
    fwrite(Vzy, 8, NumCoordinates, fp_Results);
    fwrite(Vzz, 8, NumCoordinates, fp_Results);

    fclose(fp_Results);

    return 0;
}
