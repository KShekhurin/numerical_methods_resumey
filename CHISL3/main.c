#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct slae {
    int size;
    double** A;
    double* B;
} SLAE_t;

void multMxV(int size, double** A, double* X, double* newX) {
    for (int i = 0; i < size; i++) {
        newX[i] = 0;
        for (int j = 0; j < size; j++) {
            newX[i] += A[i][j] * X[j];
        }
    }
}

void subVec(int size, double* A, double* B, double* new) {
    for (int i = 0; i < size; i++) {
        new[i] = A[i] - B[i];
    }
}

void numberMultVec(int size, double num, double* A, double* new) {
    for (int i = 0; i < size; i++) {
        new[i] = num * A[i];
    }
}

double dotProd(int size, double* A, double* B) {
    double res = 0;
    for (int i = 0; i < size; i++) {
        res += A[i] * B[i];
    }

    return res;
}

void printVec(int size, double* vec) {
    for (int i = 0; i < size; i++) {
        printf("%lf ", vec[i]);
    }
    printf("\n");
}

void logVector(FILE* f, int size, double* vec, int iter) {
    fprintf(f, "%i;", iter);
    for(int i = 0; i < size; i++) {
        fprintf(f, "%.16lf;", vec[i]);
    }
    fprintf(f, "\n");
}

double gradEndsCalc(double a_k, SLAE_t eq, double* A_g_k, double* g_k) {
    return fabs(a_k) * sqrt(dotProd(eq.size, A_g_k, g_k));
}

int gradSol(SLAE_t eq, double condA, double eps, double* X_0,
            FILE* logFile, char to_log) {
    double M = (condA - 1) / (condA + 1);
    M = M / (1 - M);
    double *g_k = malloc(sizeof(double) * eq.size);
    double *A_g_k = malloc(sizeof(double) * eq.size);
    int iter = 1;

    multMxV(eq.size, eq.A, X_0, g_k);
    for(int i = 0; i < eq.size; i++) {
        g_k[i] = 2*(g_k[i] - eq.B[i]);
    }
    multMxV(eq.size, eq.A, g_k, A_g_k);
    double a_k = dotProd(eq.size, g_k, g_k) /
            (2.0 * dotProd(eq.size, A_g_k, g_k));

    while(gradEndsCalc(a_k, eq, A_g_k, g_k) >= eps / M) {
        //printf("%0.16lf\n", a_k);

        for(int i = 0; i < eq.size; i++) {
            X_0[i] -= g_k[i] * a_k;
        }
        //printVec(eq.size, g_k);
        //printVec(eq.size, X_0);
        if(to_log)
            logVector(logFile, eq.size, X_0, iter);

        multMxV(eq.size, eq.A, X_0, g_k);
        for(int i = 0; i < eq.size; i++) {
            g_k[i] = 2*(g_k[i] - eq.B[i]);
        }
        multMxV(eq.size, eq.A, g_k, A_g_k);

        a_k = dotProd(eq.size, g_k, g_k) /
            2.0 / dotProd(eq.size, A_g_k, g_k);
        iter++;
    }

    free(g_k);
    free(A_g_k);

    return 1;
}

void Sol(char* filename, char* resName, double eps) {
    FILE* f = fopen(filename, "rb");
    FILE* r = fopen(resName, "wb");
    int32_t n;
    double cond;

    fread(&n, sizeof(int32_t), 1, f);
    fread(&cond, sizeof(double), 1, f);
    double** A = malloc(n * sizeof(double*));
    double* B = malloc(n * sizeof(double));
    double* sol = calloc(n, sizeof(double));

    for(int i = 0; i < n; i++) {
        A[i] = malloc(n * sizeof(double));
        for(int j = 0; j < n; j++) {
            fread(&A[i][j], sizeof(double), 1, f);
        }
    }

    for(int i = 0; i < n; i++) {
        fread(&B[i], sizeof(double), 1, f);
    }

    fclose(f);

    SLAE_t eq = {n, A, B};

    gradSol(eq, cond, eps, sol, r, 1);

    for(int i = 0; i < n; i++) {
        printf("%lf\n", sol[i]);
    }
    printf("\n");

    fclose(r);
    free(A);
    free(B);
    free(sol);
}

void testConvergence(char* in, char* out) {
    FILE* f = fopen(in, "rb");
    FILE* r = fopen(out, "w");
    int32_t n;
    double cond;

    fread(&n, sizeof(int32_t), 1, f);
    fread(&cond, sizeof(double), 1, f);
    double** A = malloc(n * sizeof(double*));
    double* B = malloc(n * sizeof(double));
    double* sol = calloc(n, sizeof(double));

    for(int i = 0; i < n; i++) {
        A[i] = malloc(n * sizeof(double));
        for(int j = 0; j < n; j++) {
            fread(&A[i][j], sizeof(double), 1, f);
        }
    }

    for(int i = 0; i < n; i++) {
        fread(&B[i], sizeof(double), 1, f);
    }

    fclose(f);

    SLAE_t eq = {n, A, B};

    for(int i = 0; i <= 10; i++) {
        for(int j = 0; j < n; j++) {
            sol[j] = 0;
        }
        double eps = pow(10, -1*i);

        gradSol(eq, cond, eps, sol, r, 0);

        fprintf(r, "%.16lf;", eps);
        for(int i = 0; i < n; i++) {
            fprintf(r, "%.16lf;", sol[i]);
        }
        fprintf(r, "\n");
    }

    fclose(r);
    free(A);
    free(B);
    free(sol);
}


int main() {
    Sol("grad10_SRC.bin", "grad10_Err_RES.csv", 1e-12);
    Sol("grad3_SRC.bin", "grad3_Err_RES.csv", 1e-12);
    Sol("grad15_SRC.bin", "grad15_Err_RES.csv", 1e-12);
    Sol("grad_VIS_SRC.bin", "grad_VIS_RES.csv", 1e-5);
    //Sol("grad_TEST.bin", "grad_TEST.csv", 1e-15);

    testConvergence("grad3_SRC.bin", "grad3_Conv_RES.csv");
    testConvergence("grad10_SRC.bin", "grad10_Conv_RES.csv");
    testConvergence("grad15_SRC.bin", "grad15_Conv_RES.csv");
    return 0;
}