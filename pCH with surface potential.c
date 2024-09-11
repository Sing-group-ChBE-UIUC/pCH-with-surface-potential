#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
pid_t getpid(void);

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define NDIV (1 + IMM1 / NTAB)
#define SWAP(a, b) \
    tempr = (a);   \
    (a) = (b);     \
    (b) = tempr

void fourn(double data[], unsigned long nn[], int ndim, int isign);
float ran1(long *idum);
float gasdev(long *idum);
long initRan();
void FirstOrderApprox(int SpeciesType, int CurrentMaxN, int N0, int Z, double Chi, double **PhiX, double **PhiY, double *SumPhi_SameType, double *SumPhi_DiffType, double **sgrad, double **Mu, double **lap, double **gradgrad, double dt, double Kappa, double Alpha, double **Phi_temp, double *SumPhi_temp, double *PhiL, double V, double Gamma);
void phi_normalization(int Z, double *Phi, double PhiA_i);
void PolymerInitialization(double **Phi_A, double **Phi_B, double *SumPhiA, double *SumPhiB, int Z, long *idum, double PhiA_i, double *PhiL);
void SaveInput(int Class, double Gamma, double V, double dt, double Chi, int Max_NA_sys, int Max_NB_sys, double PhiA_i, double PhiB_i, double k_rxn_A, double k_rxn_B, double Kappa, int N0, double Alpha, int t_max, int Z, int t_save, long *idum);
void Save(int Class, double *quantity, double *TotalPhi, int N, int Z, int t, char Content[]);
int SaveOutput(double *SumPhiA, double *SumPhiB, double **Phi_A, double **Phi_B, int Class, int t, int Max_NA_current, int Max_NB_current, int Z);
void LinearStepGrowthPolymerization(double **Phi, int Z2, double k_rxn, int Max_N_current);

int main(int argc, const char *argv[])
// int main()
{
    int i, k, n, j;
    int NA_max; // Max degree of polymerization possible in the system
    int NB_max; // Max degree of polymerization possible in the system
    int Z;      // number of data points on one side

    double sumphiZ;

    // random seed
    long *idum = malloc(sizeof(long));
    *idum = initRan();
    // *idum = 1;

    int t;
    int NA_max_cap = 1;
    int NB_max_cap = 1;
    int Class;
    int t_max; // total iterations
    double Kappa, Chi, Alpha;
    double k_A, k_B;
    double Kappatemp;
    double dt;
    double PhiA_i, PhiB_i;
    int t_save;
    int N0;
    double V; // Potential at wall (will decay away from the wall)
    double Gamma; // Potential decay strength

    char *SaveStr1 = malloc(sizeof(char) * 100); // memory for output file name
    FILE *Save1;

    sscanf(argv[1], "%d", &Class); // read from command line argument. You will execute your code as './a.out %d', where the integer is the trajectory #

    // Class = 1;
    // sprintf(SaveStr1, "%d-Input.txt", Class);
    sprintf(SaveStr1, "Input.txt");
    Save1 = fopen(SaveStr1, "r");

    fscanf(Save1, "dt = %lf\n", &dt);
    fscanf(Save1, "Chi = %lf\n", &Chi);
    fscanf(Save1, "NA_max = %d\n", &NA_max);
    fscanf(Save1, "NB_max = %d\n", &NB_max);
    fscanf(Save1, "PhiA_i = %lf\n", &PhiA_i);
    fscanf(Save1, "PhiB_i = %lf\n", &PhiB_i);
    fscanf(Save1, "k_rxn_A = %lf\n", &k_A);
    fscanf(Save1, "k_rxn_B = %lf\n", &k_B);
    fscanf(Save1, "V = %lf\n", &V);
    fscanf(Save1, "Gamma = %lf\n", &Gamma);
    fscanf(Save1, "Kappa = %lf\n", &Kappa);
    fscanf(Save1, "N_0 = %d\n", &N0);
    fscanf(Save1, "Alpha = %lf\n", &Alpha);
    fscanf(Save1, "Z = %d\n", &Z);
    fscanf(Save1, "t_max = %d\n", &t_max);
    fscanf(Save1, "t_save = %d\n", &t_save);
    fclose(Save1);

    int Z2 = Z * Z;
    int check;

    int t_ramp = 10000;

    double double_temp, TotalSum;

    // initialization for numerical calculation
    double **Phi_A, **Phi_B;           // phi A and phi B in real space ( format -> Phi_A[DoP][Grid Index] )
    double **Phi_A_temp, **Phi_B_temp; // phi A and phi B in real space ( format -> Phi_A[DoP][Grid Index] )
    double **sgradA, **muA, **sgradB, **muB, **lapA, **lapB, **gradgradA, **gradgradB;

    Phi_A = (double **)calloc(NA_max, sizeof(double *));
    Phi_A_temp = (double **)calloc(NA_max, sizeof(double *));
    sgradA = (double **)calloc(NA_max, sizeof(double *));
    muA = (double **)calloc(NA_max, sizeof(double *));
    lapA = (double **)calloc(NA_max, sizeof(double *));
    gradgradA = (double **)calloc(NA_max, sizeof(double *));

    for (i = 0; i < NA_max; i++)
    {
        Phi_A[i] = (double *)calloc(Z2, sizeof(double));
        Phi_A_temp[i] = (double *)calloc(Z2, sizeof(double));
        sgradA[i] = (double *)calloc(Z2, sizeof(double));
        muA[i] = (double *)calloc(Z2, sizeof(double));
        lapA[i] = (double *)calloc(Z2, sizeof(double));
        gradgradA[i] = (double *)calloc(Z2, sizeof(double));
    }

    Phi_B = (double **)calloc(NB_max, sizeof(double *));
    Phi_B_temp = (double **)calloc(NB_max, sizeof(double *));
    sgradB = (double **)calloc(NB_max, sizeof(double *));
    muB = (double **)calloc(NB_max, sizeof(double *));
    lapB = (double **)calloc(NB_max, sizeof(double *));
    gradgradB = (double **)calloc(NB_max, sizeof(double *));

    for (i = 0; i < NB_max; i++)
    {
        Phi_B[i] = (double *)calloc(Z2, sizeof(double));
        Phi_B_temp[i] = (double *)calloc(Z2, sizeof(double));
        sgradB[i] = (double *)calloc(Z2, sizeof(double));
        muB[i] = (double *)calloc(Z2, sizeof(double));
        lapB[i] = (double *)calloc(Z2, sizeof(double));
        gradgradB[i] = (double *)calloc(Z2, sizeof(double));
    }

    double *SumPhiA = calloc(Z2, sizeof(double)); // Sum of Phi_A
    double *SumPhiB = calloc(Z2, sizeof(double)); // Sum of Phi_B

    double *SumPhiA_temp = calloc(Z2, sizeof(double)); // Sum of Phi_A
    double *SumPhiB_temp = calloc(Z2, sizeof(double)); // Sum of Phi_B

    double *PhiL = calloc(Z2, sizeof(double)); // Sum of Phi_A

    // This below sets the interface by representing a total phi(z), where z is the direction normal to the surface, and phi(z)<<1 as we approach the surface.
    for (i = 0; i < Z; ++i)
    {
        for (j = 0; j < Z; ++j)
        {
            PhiL[Z * i + j] = 0.25 * (1.0 + tanh(i - 4)) * (1.0 + tanh(Z - 5 - i));
            // printf("%lf ", PhiL[Z*i+j]);
        }
        // printf("\n");
    }

    for (k = 0; k <= 1; k++)
    {
        // reset all memory
        for (i = 0; i < NA_max; i++)
        {
            memset(Phi_A[i], 0, Z2 * sizeof(double));
            memset(sgradA[i], 0, Z2 * sizeof(double));
            memset(muA[i], 0, Z2 * sizeof(double));
            memset(lapA[i], 0, Z2 * sizeof(double));
            memset(gradgradA[i], 0, Z2 * sizeof(double));
        }
        memset(*Phi_A, 0, NA_max * sizeof(double *));
        memset(*sgradA, 0, NA_max * sizeof(double *));
        memset(*muA, 0, NA_max * sizeof(double *));
        memset(*lapA, 0, NA_max * sizeof(double *));
        memset(*gradgradA, 0, NA_max * sizeof(double *));

        for (i = 0; i < NB_max; i++)
        {
            memset(Phi_B[i], 0, Z2 * sizeof(double));
            memset(sgradB[i], 0, Z2 * sizeof(double));
            memset(muB[i], 0, Z2 * sizeof(double));
            memset(lapB[i], 0, Z2 * sizeof(double));
            memset(gradgradB[i], 0, Z2 * sizeof(double));
        }
        memset(*Phi_B, 0, NB_max * sizeof(double *));
        memset(*sgradB, 0, NB_max * sizeof(double *));
        memset(*muB, 0, NB_max * sizeof(double *));
        memset(*lapB, 0, NB_max * sizeof(double *));
        memset(*gradgradB, 0, NB_max * sizeof(double *));

        memset(SumPhiA, 0, Z2 * sizeof(double));
        memset(SumPhiB, 0, Z2 * sizeof(double));
        memset(SumPhiA_temp, 0, Z2 * sizeof(double));
        memset(SumPhiB_temp, 0, Z2 * sizeof(double));

        NA_max_cap = 1;
        NB_max_cap = 1;

        Class += k;
        SaveInput(Class, Gamma, V, dt, Chi, NA_max, NB_max, PhiA_i, PhiB_i, k_A, k_B, Kappa, N0, Alpha, t_max, Z, t_save, idum);
        PolymerInitialization(Phi_A, Phi_B, SumPhiA, SumPhiB, Z, idum, PhiA_i, PhiL);

        check = 0;
        // Enter Updating Loop
        for (t = 0; t <= t_max; t++)
        {

            if (t % t_save == 0)
            {
                if (SaveOutput(SumPhiA, SumPhiB, Phi_A, Phi_B, Class, t, NA_max, NB_max, Z) == 1)
                {
                    break; // move on to next simulation
                }
            }

            // Need to ramp up kappa (this is relatively quick), so we don't break the first iteration or two.
            if (t < t_ramp)
            {
                Kappatemp = Kappa * (0.1 + 0.9 * t / t_ramp);
            }
            else
            {
                Kappatemp = Kappa;
            }

            // A species
            FirstOrderApprox(0, NA_max_cap, N0, Z, Chi, Phi_A, Phi_B, SumPhiA, SumPhiB, sgradA, muA, lapA, gradgradA, dt, Kappatemp, Alpha, Phi_A_temp, SumPhiA_temp, PhiL, V, Gamma);

            // B species
            FirstOrderApprox(1, NB_max_cap, N0, Z, Chi, Phi_B, Phi_A, SumPhiB, SumPhiA, sgradB, muB, lapB, gradgradB, dt, Kappatemp, Alpha, Phi_B_temp, SumPhiB_temp, PhiL, V, Gamma);

            // rxn doubles DOP in the system until it reaches Max DOP
            if (NA_max_cap < NA_max)
            {
                NA_max_cap *= 2;
            }

            // rxn doubles DOP in the system until it reaches Max DOP
            if (NB_max_cap < NB_max)
            {
                NB_max_cap *= 2;
            }

            // Reaction in all grids
            LinearStepGrowthPolymerization(Phi_A, Z2, k_A, NA_max_cap);
            LinearStepGrowthPolymerization(Phi_B, Z2, k_B, NB_max_cap);

            memset(SumPhiA, 0, Z2 * sizeof(double));
            memset(SumPhiB, 0, Z2 * sizeof(double));

            for (i = 0; i < Z2; ++i)
            {
                double_temp = 0.0;
                for (n = 0; n < NA_max_cap; ++n)
                {
                    SumPhiA[i] += Phi_A[n][i];
                }
                for (n = 0; n < NB_max_cap; ++n)
                {
                    SumPhiB[i] += Phi_B[n][i];
                }
                double_temp = SumPhiA[i] + SumPhiB[i];

                if (isnan(double_temp))
                {
                    check = SaveOutput(SumPhiA, SumPhiB, Phi_A, Phi_B, Class, t, NA_max, NB_max, Z);
                    check = 1;
                    break; // end i loop
                }
            }

            if (check == 1)
            {
                break; // end t loop
            }
        }
    }
    return 0;
}

int SaveOutput(double *SumPhiA, double *SumPhiB, double **Phi_A, double **Phi_B, int Class, int t, int Max_NA_current, int Max_NB_current, int Z)
{
    int i, j, N;
    int Z2 = Z * Z;
    double *TotalPhi = calloc(Z2, sizeof(double)); // Sum of Phi_A

    char *SaveStr2 = malloc(sizeof(char) * 100); // memory for output file name
    FILE *Save2;

    sprintf(SaveStr2, "%d-Psi-Overall.txt", Class);
    Save2 = fopen(SaveStr2, "a"); // open output file (overall Psi)

    // initialize variables
    double MixingParam = 0; // Degree of mixing
    double AvgPhiA = 0;     // Average PhiA
    double AvgPhiB = 0;     // Average PhiB
    double Psi = 0;         // Degree of phase separation
    double sum = 0;
    double Weight2 = 0;
    double Weight4 = 0;
    double Weight10 = 0;

    // calculate variables
    for (i = 0; i < Z2; i++) // loop over all grid space
    {
        TotalPhi[i] = SumPhiA[i] + SumPhiB[i];
        Weight2 = TotalPhi[i] * TotalPhi[i];
        Weight4 = Weight2 * Weight2;
        Weight10 = Weight4 * Weight4 * Weight2;
        
        MixingParam += (SumPhiA[i] * SumPhiB[i]) * Weight10;
        sum += Weight10;
        AvgPhiA += SumPhiA[i] * Weight10;
        AvgPhiB += SumPhiB[i] * Weight10;
    }

    // AvgPhiA /= (double)Z2;                         // average over space
    // AvgPhiB /= (double)Z2;                         // average over space
    MixingParam *= sum / (AvgPhiA * AvgPhiB); // normalize and average over space
    Psi = 1.0 - MixingParam;                  // calculate Psi

    fprintf(Save2, "%d %lf\n", t, Psi); // write output file (overall Psi)
    fclose(Save2);                      // close output file (overall Psi)

    // generate output file for Phi A over different N
    for (N = 0; N < Max_NA_current; N++)
    {
        Save(Class, Phi_A[N], TotalPhi, N, Z, t, "PhiA");
    }

    // // generate output file for Phi B over different N
    // for (N = 0; N < Max_NB_current; N++)
    // {
    //     Save(Class, Phi_B[N], TotalPhi, N, Z, t, "PhiB");
    // }

    Save(Class, SumPhiA, TotalPhi, -1, Z, t, "SumPhiA");   // generate output file for sum of Phi A
    Save(Class, TotalPhi, TotalPhi, -1, Z, t, "TotalPhi"); // generate output file for Total Phi
    // Save(Class, SumPhiB, TotalPhi, -1, Z, t, "SumPhiB"); // generate output file for sum of Phi B

    free(TotalPhi);

    if (Psi > 0.45)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void Save(int Class, double *quantity, double *TotalPhi, int N, int Z, int t, char Content[])
{
    int i, j, ij;

    // memory for file name
    char *FileName = malloc(sizeof(char) * 100);
    FILE *save;

    sprintf(FileName, "%d-%s-N%d-t%d.txt", Class, Content, N + 1, t); // set file name
    save = fopen(FileName, "w");                                      // open file

    // Over all grid space
    for (i = 0; i < Z; i++)
    {
        for (j = 0; j < Z; j++)
        {
            ij = j * Z + i; // spatial location (or index)
            // fprintf(save, "%lf ", quantity[ij] / TotalPhi[ij]); // write file
            fprintf(save, "%lf ", quantity[ij]); // write file
        }
        fprintf(save, "\n"); // switch row in written file
    }
    fclose(save); // close file

    free(FileName); // free memory
}

// Input save function
void SaveInput(int Class, double Gamma, double V, double dt, double Chi, int Max_NA_sys, int Max_NB_sys, double PhiA_i, double PhiB_i, double k_rxn_A, double k_rxn_B, double Kappa, int N0, double Alpha, int t_max, int Z, int t_save, long *idum)
{
    double C = sqrt(1.0 / (Kappa * N0 * 6.0 * PhiA_i * PhiB_i));
    double L0 = C * N0; // grid distance units of Rg

    // memory for file name
    char *FileName = malloc(sizeof(char) * 100);
    FILE *save;

    sprintf(FileName, "%d-Input.txt", Class); // set file name
    save = fopen(FileName, "w");              // open file

    // save variables
    fprintf(save, "dt = %lf\n", dt);
    fprintf(save, "Chi = %lf\n", Chi);
    fprintf(save, "NA_max = %d\n", Max_NA_sys);
    fprintf(save, "NB_max = %d\n", Max_NB_sys);
    fprintf(save, "PhiA_i = %lf\n", PhiA_i);
    fprintf(save, "PhiB_i = %lf\n", PhiB_i);
    fprintf(save, "k_rxn_A = %lf\n", k_rxn_A);
    fprintf(save, "k_rxn_B = %lf\n", k_rxn_B);
    fprintf(save, "V = %lf\n", V);
    fprintf(save, "Gamma = %lf\n", Gamma);
    fprintf(save, "N_0 = %d\n", N0);
    fprintf(save, "Kappa = %lf\n", Kappa);
    fprintf(save, "L0 = %lf\n", L0);
    fprintf(save, "Alpha = %lf\n", Alpha);
    fprintf(save, "Z = %d\n", Z);
    fprintf(save, "t_max = %d\n", t_max);
    fprintf(save, "t_save = %d\n", t_save);
    fprintf(save, "\nseed = %ld\n", *idum);

    fclose(save); // close file

    free(FileName); // free memory
}

void FirstOrderApprox(int SpeciesType, int CurrentMaxN, int N0, int Z, double Chi, double **PhiX, double **PhiY, double *SumPhi_SameType, double *SumPhi_DiffType, double **sgrad, double **Mu, double **lap, double **gradgrad, double dt, double Kappa, double Alpha, double **Phi_temp, double *SumPhi_temp, double *PhiL, double V, double Gamma)
{
    int i, j, g0, gx1, gxm1, gy1, gym1, gx2, gy2, gxm2, gym2, N;
    int gxm3, gx3;
    double dphix, dphiy, dmux, dmuy;
    double V0;
    double PotentialExp = Gamma;
    double WallPotential;
    double phiconst;

    if (SpeciesType == 0) // A type
    {
        V0 = 0.0;
    }
    else // B type
    {
        V0 = V;
    }

    for (i = 0; i < Z; i++)
    {
        if (i < 3)
        {
            WallPotential = 0;
        }
        else
        {
            WallPotential = V0 * exp(-PotentialExp * (i - 3.0));
        }

        for (j = 0; j < Z; j++)
        {
            g0 = Z * i + j;
            gx1 = Z * ((i + 1) % Z) + j;
            gx2 = Z * ((i + 2) % Z) + j;
            gxm1 = Z * ((i - 1 + Z) % Z) + j;
            gxm2 = Z * ((i - 2 + Z) % Z) + j;
            gy1 = Z * i + ((j + 1) % Z);
            gy2 = Z * i + ((j + 2) % Z);
            gym1 = Z * i + ((j - 1 + Z) % Z);
            gym2 = Z * i + ((j - 2 + Z) % Z);

            for (N = 0; N < CurrentMaxN; N++)
            {
                // Keep it simple! If I change the constraint to be at a total phi(z) that is <1 near the surfaces, then I don't need to modify the Cahn-Hilliard update.

                // x-dir
                sgrad[N][g0] = (16.0 * (PhiX[N][gx1] / (SumPhi_SameType[gx1] + SumPhi_DiffType[gx1]) + PhiX[N][gxm1] / (SumPhi_SameType[gxm1] + SumPhi_DiffType[gxm1])) - (PhiX[N][gx2] / (SumPhi_SameType[gx2] + SumPhi_DiffType[gx2]) + PhiX[N][gxm2] / (SumPhi_SameType[gxm2] + SumPhi_DiffType[gxm2])) - 30.0 * PhiX[N][g0] / (SumPhi_SameType[g0] + SumPhi_DiffType[g0])) / 12.0;

                // y-dir (additive)
                sgrad[N][g0] += (16.0 * (PhiX[N][gy1] / (SumPhi_SameType[gy1] + SumPhi_DiffType[gy1]) + PhiX[N][gym1] / (SumPhi_SameType[gym1] + SumPhi_DiffType[gym1])) - (PhiX[N][gy2] / (SumPhi_SameType[gy2] + SumPhi_DiffType[gy2]) + PhiX[N][gym2] / (SumPhi_SameType[gym2] + SumPhi_DiffType[gym2])) - 30.0 * PhiX[N][g0] / (SumPhi_SameType[g0] + SumPhi_DiffType[g0])) / 12.0;
            }

            // HERE is where I include the phi(z) constraint (as PhiL[g0] in the below expression).
            for (N = 0; N < CurrentMaxN; ++N)
            {
                phiconst = Alpha * (SumPhi_SameType[g0] + SumPhi_DiffType[g0] - PhiL[g0]);
                // if (phiconst > 2.0)
                // {
                //     phiconst = 2.0;
                // }
                // if (phiconst < -2.0)
                // {
                //     phiconst = -2.0;
                // }
                // Mu[N][g0] += phiconst;
                Mu[N][g0] = (log(PhiX[N][g0]) + 1.0) / ((double)N + 1.0) + Chi * SumPhi_DiffType[g0] + phiconst - Kappa * sgrad[N][g0] / 2.0 + WallPotential;

                Phi_temp[N][g0] = PhiX[N][g0];
            }
            SumPhi_temp[g0] = SumPhi_SameType[g0];
        }
    }

    for (i = 0; i < Z; i++)
    {
        for (j = 0; j < Z; j++)
        {
            g0 = Z * i + j;
            gx1 = Z * ((i + 1) % Z) + j;
            gx2 = Z * ((i + 2) % Z) + j;
            gxm1 = Z * ((i - 1 + Z) % Z) + j;
            gxm2 = Z * ((i - 2 + Z) % Z) + j;
            gy1 = Z * i + ((j + 1) % Z);
            gy2 = Z * i + ((j + 2) % Z);
            gym1 = Z * i + ((j - 1 + Z) % Z);
            gym2 = Z * i + ((j - 2 + Z) % Z);

            for (N = 0; N < CurrentMaxN; ++N)
            {
                // x-dir
                lap[N][g0] = (16.0 * (Mu[N][gx1] + Mu[N][gxm1]) - (Mu[N][gx2] + Mu[N][gxm2]) - 30.0 * Mu[N][g0]) / 12.0;
                dphix = (PhiX[N][gx1] - PhiX[N][gxm1]) / 2.0;
                dmux = (Mu[N][gx1] - Mu[N][gxm1]) / 2.0;

                // y-dir (additive)
                lap[N][g0] += (16.0 * (Mu[N][gy1] + Mu[N][gym1]) - (Mu[N][gy2] + Mu[N][gym2]) - 30.0 * Mu[N][g0]) / 12.0;

                dphiy = (PhiX[N][gy1] - PhiX[N][gym1]) / 2.0;
                dmuy = (Mu[N][gy1] - Mu[N][gym1]) / 2.0;
                gradgrad[N][g0] = (dphix * dmux + dphiy * dmuy) / PhiX[N][g0];
                if (gradgrad[N][g0] > 1.0)
                    gradgrad[N][g0] = 1.0; // Might not need this - will depend on kappa, which needs to be sufficiently big
                if (gradgrad[N][g0] < -1.0)
                    gradgrad[N][g0] = -1.0;
            }
        }
    }

    for (i = 0; i < Z; ++i)
    {
        for (j = 0; j < Z; ++j)
        {
            g0 = Z * i + j;
            for (N = 0; N < CurrentMaxN; ++N)
            {
                PhiX[N][g0] = log(PhiX[N][g0]);
                PhiX[N][g0] += N0 * dt * (lap[N][g0] + gradgrad[N][g0]) / ((double)N + 1.0);
                PhiX[N][g0] = exp(PhiX[N][g0]);
            }
        }
    }
}

void LinearStepGrowthPolymerization(double **Phi, int Z2, double k_rxn, int Max_N_current)
{
    int i, j, k;
    double Generation, Consumption;

    // Temp variable
    double *delta_Phi = calloc(Max_N_current, sizeof(double));

    for (k = 0; k < Z2; k++)
    {
        // Upto Current max DOP
        for (i = 0; i < Max_N_current; i++)
        {
            Generation = 0.0;

            // Generation
            for (j = 0; j < i; j++)
            {
                // N = i-1 chains formed based on concentration ~ phi/DOP
                Generation += Phi[j][k] / (j + 1.0) * Phi[i - j - 1][k] / (i - j);
            }

            Consumption = 0;

            // i-mer is consumed with all possible j-mers in the system except for N = MaxN chain since it is not consumed.
            for (j = 0; j < Max_N_current - i - 1; j++)
            {
                // Concentration ~ phi/DOP
                Consumption += Phi[j][k] / (j + 1.0);
            }

            // Convert concentration to phi && calc. delta phi
            delta_Phi[i] = k_rxn * (0.5 * Generation * (i + 1.0) - Phi[i][k] * Consumption);
        }

        // apply changes
        for (i = 0; i < Max_N_current; i++)
        {
            // Apply change
            Phi[i][k] += delta_Phi[i];
        }
    }

    free(delta_Phi);
}

void PolymerInitialization(double **Phi_A, double **Phi_B, double *SumPhiA, double *SumPhiB, int Z, long *idum, double PhiA_i, double *PhiL)
{
    int i, j;
    double factor, factor1;

    factor = 0.0;
    factor1 = 0.0;
    for (i = 0; i < Z * Z; ++i)
    {
        factor1 += PhiL[i];
    }
    factor1 /= (double)(Z * Z); /**/

    // Monomer A initialization
    for (i = 0; i < Z; ++i)
    {
        for (j = 0; j < Z; j++)
        {
            Phi_A[0][i * Z + j] = (PhiA_i + 0.001 * (double)ran1(idum)) * PhiL[i * Z + j]; //(PhiA_i + PhiA_i * 0.001 * (double)ran1(idum))*factor;
            factor += Phi_A[0][i * Z + j];
        }
    }
    factor /= (double)(Z * Z);
    // printf("%lf\n", factor);
    factor /= (PhiA_i * factor1);
    // printf("%lf\n", factor);

    // Monomer A normalization such that <PhiA> = PhiA_i
    // phi_normalization(Z, Phi_A[0], PhiA_i);

    for (i = 0; i < Z * Z; ++i)
    {
        Phi_A[0][i] /= factor;
        Phi_B[0][i] = PhiL[i] - Phi_A[0][i]; // Monomer B
        SumPhiA[i] = Phi_A[0][i];            // Sum of PhiA
        SumPhiB[i] = Phi_B[0][i];            // Sum of PhiB
    }
}

// Phi normalization function (only used once at initialization)
void phi_normalization(int Z, double *Phi, double PhiA_i)
{
    double Sum = 0;
    int i, j, ij;

    // get sum of PhiA
    for (i = 0; i < Z; i++)
    {
        for (j = 0; j < Z; j++)
        {
            ij = Z * j + i;
            Sum += Phi[ij];
        }
    }
    Sum /= (Z * Z); // average sum of PhiA

    // Normalize PhiA such that <PhiA> = PhiA_i
    for (i = 0; i < Z; i++)
    {
        for (j = 0; j < Z; j++)
        {
            ij = Z * j + i;

            Phi[ij] *= PhiA_i / Sum;
        }
    }
}

void fourn(double data[], unsigned long nn[], int ndim, int isign)
// Replaces data by its ndim-dimensional discrete Fourier transform, if isign is input as 1.
// nn[1..ndim] is an integer array containing the lengths of each dimension (number of complex
// values), which MUST all be powers of 2. data is a real array of length twice the product of
// these lengths, in which the data are stored as in a multidimensional complex array: real and
// imaginary parts of each element are in consecutive locations, and the rightmost index of the
// array increases most rapidly as one proceeds along data. For a two-dimensional array, this is
// equivalent to storing the array by rows. If isign is input as ?��1, data is replaced by its inverse
// transform times the product of the lengths of all dimensions.
{
    int idim;
    unsigned long i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
    unsigned long ibit, k1, k2, n, nprev, nrem, ntot;
    float tempi, tempr;
    double theta, wi, wpi, wpr, wr, wtemp;         // Double precision for trigonometric recurrences.
    for (ntot = 1, idim = 1; idim <= ndim; idim++) // Compute total number of complex values.
        ntot *= nn[idim];
    nprev = 1;
    for (idim = ndim; idim >= 1; idim--)
    { // Main loop over the dimensions.
        n = nn[idim];
        nrem = ntot / (n * nprev);
        ip1 = nprev << 1;
        ip2 = ip1 * n;
        ip3 = ip2 * nrem;
        i2rev = 1;
        for (i2 = 1; i2 <= ip2; i2 += ip1)
        { // This is the bit-reversal section of the routine.
            if (i2 < i2rev)
            {
                for (i1 = i2; i1 <= i2 + ip1 - 2; i1 += 2)
                {
                    for (i3 = i1; i3 <= ip3; i3 += ip2)
                    {
                        i3rev = i2rev + i3 - i2;
                        SWAP(data[i3], data[i3rev]);
                        SWAP(data[i3 + 1], data[i3rev + 1]);
                    }
                }
            }
            ibit = ip2 >> 1;
            while (ibit >= ip1 && i2rev > ibit)
            {
                i2rev -= ibit;
                ibit >>= 1;
            }
            i2rev += ibit;
        }
        ifp1 = ip1; // Here begins the Danielson-Lanczos section of the routine.
        while (ifp1 < ip2)
        {
            ifp2 = ifp1 << 1;
            theta = isign * 6.28318530717959 / (ifp2 / ip1); // Initialize for the trig. recurrence.
            wtemp = sin(0.5 * theta);
            wpr = -2.0 * wtemp * wtemp;
            wpi = sin(theta);
            wr = 1.0;
            wi = 0.0;
            for (i3 = 1; i3 <= ifp1; i3 += ip1)
            {
                for (i1 = i3; i1 <= i3 + ip1 - 2; i1 += 2)
                {
                    for (i2 = i1; i2 <= ip3; i2 += ifp2)
                    {
                        k1 = i2; // Danielson-Lanczos formula:
                        k2 = k1 + ifp1;
                        tempr = (float)wr * data[k2] - (float)wi * data[k2 + 1];
                        tempi = (float)wr * data[k2 + 1] + (float)wi * data[k2];
                        data[k2] = data[k1] - tempr;
                        data[k2 + 1] = data[k1 + 1] - tempi;
                        data[k1] += tempr;
                        data[k1 + 1] += tempi;
                    }
                }
                wr = (wtemp = wr) * wpr - wi * wpi + wr; // Trigonometric recurrence.
                wi = wi * wpr + wtemp * wpi + wi;
            }
            ifp1 = ifp2;
        }
        nprev *= n;
    }
}

/// Numerical recipes for random number generation - this one is between 0 and 1
float ran1(long *idum)
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0)
    {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; --j)
        {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0)
        *idum += IM1;
    k = idum2 / IQ2;
    if (*idum < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

/// Numerical recipes for random number generation - this one is Gaussian distributed around 0, with a variance of 1.
float gasdev(long *idum)
{
    float ran1(long *idum);
    static int iset = 0;
    static float gset;
    float fac, rsq, v1, v2;

    if (*idum < 0)
        iset = 0;
    if (iset == 0)
    {
        do
        {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}

// random seed initialization
long initRan()
{
    // time_t seconds;
    // time(&seconds);
    // return -1*(unsigned long)(seconds/12345); This is bad.  :disappointed:

    // This will hopefully allow us to have a unique seed even if executed multiple times a second-Got from Mike
    // http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a = a - b;
    a = a - c;
    a = a ^ (c >> 13);
    b = b - c;
    b = b - a;
    b = b ^ (a << 8);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 13);
    a = a - b;
    a = a - c;
    a = a ^ (c >> 12);
    b = b - c;
    b = b - a;
    b = b ^ (a << 16);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 5);
    a = a - b;
    a = a - c;
    a = a ^ (c >> 3);
    b = b - c;
    b = b - a;
    b = b ^ (a << 10);
    c = c - a;
    c = c - b;
    c = c ^ (b >> 15);
    return c % 1000000000; // careful here.  Another 0 might break the ran1 (long long instead of just long)
}
