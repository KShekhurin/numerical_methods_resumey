#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct matrix {
    int size;
    double** A;
    double* B;
} Matrix_t;

double** initMat(int n) {
    double** A = calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++) {
        A[i] = calloc(n, sizeof(double));
    }

    return A;
}

double dotProduct(int n, double* A, double* B) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += A[i] * B[i];
    }
    return result;
}

double norm_2(int n, double* A) {
    return sqrt(dotProduct(n, A, A));
}

double norm_INF(int n, double* V) {
    double res = 0;
    for (int i = 0; i < n; i++) {
        if(res < fabs(V[i]))
            res = fabs(V[i]);
    }
    return res;
}

void LU_decomp(int n, double** A, double** L, double** U) {
    for(int m = 0; m < n; m++) {
        L[m][m] = 1;

        for(int j = m; j < n; j++) {
            U[m][j] = A[m][j];
            for(int k = 0; k < m; k++) {
                U[m][j] -= L[m][k] * U[k][j];
            }
        }

        for(int i = m + 1; i < n; i++) {
            L[i][m] = A[i][m];
            for(int k = 0; k < m; k++) {
                L[i][m] -= L[i][k] * U[k][m];
            }
            L[i][m] /= U[m][m];
        }
    }
};

void solve_LTRIAG(int n, double** L, double* B, double* X) {
    for(int i = 0; i < n; i++) {
        X[i] = B[i];
        for(int j = 0; j < i; j++) {
            X[i] -= L[i][j] * X[j];
        }
    }
}

void solve_UTRIAG(int n, double** U, double* B, double* X) {
    for(int i = n - 1; i >= 0; i--) {
        X[i] = B[i];
        for(int j = i + 1; j < n; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
}

void scalar_mult(int n, double* V, double mu) {
    for(int i = 0; i < n; i++) {
        V[i] *= mu;
    }
}

void vlog(FILE* f, double lambda, int iter) {
    fprintf(f, "%i;%.16lf\n", iter, lambda);
}

void MxV(int size, double** A, double* X, double* out) {
    for(int i = 0; i < size; i++) {
        out[i] = 0;
        for(int j = 0; j < size; j++) {
            out[i] += A[i][j] * X[j];
        }
    }
}

double approx(Matrix_t eq, double* X, double lambda) {
    double* X1 = calloc(eq.size, sizeof(double));
    double* X2 = calloc(eq.size, sizeof(double));

    for(int i = 0; i < eq.size; i++) {
        X2[i] = X[i];
    }

    MxV(eq.size, eq.A, X, X1);
    scalar_mult(eq.size, X2, lambda);

    for(int i = 0; i < eq.size; i++) {
        X1[i] -= X2[i];
    }

    double res = norm_2(eq.size, X1) / norm_2(eq.size, X);

    free(X1);
    free(X2);

    //printf("Error: %lf\n", res);

    return res;
}

double IIM(Matrix_t eq, double* X0, double lambda, double eps, FILE* f, char to_log) {
    double** L = initMat(eq.size);
    double** U = initMat(eq.size);
    double** A = initMat(eq.size);

    double* X_prev = calloc(eq.size, sizeof(double));
    double* X_TMP = calloc(eq.size, sizeof(double));

    int iter = 0;
    double lambda_approx = lambda;

    for(int i = 0; i < eq.size; i++) {
        for(int j = 0; j < eq.size; j++) {
            A[i][j] = eq.A[i][j];
            if(i == j)
                A[i][j] -= lambda;
        }
    }
    LU_decomp(eq.size, A, L, U);

    do {
        double mu = norm_INF(eq.size, X0);

        for(int j = 0; j < eq.size; j++) {
            X_prev[j] = X0[j];
        }

        scalar_mult(eq.size, X0, 1 / mu);
        solve_LTRIAG(eq.size, L, X0, X_TMP);
        solve_UTRIAG(eq.size, U, X_TMP, X0);

        iter++;

        double shift = 0;
        for(int i = 0; i < eq.size; i++) {
            shift += X_prev[i] / X0[i];
        }
        shift /= eq.size;
        lambda_approx = lambda + shift / norm_INF(eq.size, X_prev);

        if(to_log)
            vlog(f, lambda_approx, iter);
        //printf("%lf\n", lambda_approx);
    } while(approx(eq, X0, lambda_approx) > eps);

    printf("approx: %lf\n", lambda_approx);

    free(X_prev);
    free(X_TMP);

    for(int i = 0; i < eq.size; i++) {
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);

    return lambda_approx;
}

void Sol(char* filename, char* resName, double eps) {
    FILE* f = fopen(filename, "rb");
    FILE* r = fopen(resName, "wb");
    int32_t n;
    double cond;
    double lambda;

    fread(&n, sizeof(int32_t), 1, f);
    fread(&lambda, sizeof(double), 1, f);
    double** A = malloc(n * sizeof(double*));
    double* sol = calloc(n, sizeof(double));
    for(int i = 0; i < n; i++) {
        sol[i] = 1;
    }

    for(int i = 0; i < n; i++) {
        A[i] = malloc(n * sizeof(double));
        for(int j = 0; j < n; j++) {
            fread(&A[i][j], sizeof(double), 1, f);
        }
    }

    fclose(f);

    Matrix_t eq = {n, A, NULL};

    IIM(eq, sol, lambda, eps, r, 1);

    // for(int i = 0; i < n; i++) {
    //     printf("%lf\n", sol[i]);
    // }
    // printf("\n");

    fclose(r);
    free(A);
    free(sol);
}

void Conv(char* filename, char* resName) {
    FILE* f = fopen(filename, "rb");
    FILE* r = fopen(resName, "w");
    int32_t n;
    double cond;
    double lambda;

    fread(&n, sizeof(int32_t), 1, f);
    fread(&lambda, sizeof(double), 1, f);
    double** A = malloc(n * sizeof(double*));
    double* sol = calloc(n, sizeof(double));
    for(int i = 0; i < n; i++) {
        sol[i] = 1;
    }

    for(int i = 0; i < n; i++) {
        A[i] = malloc(n * sizeof(double));
        for(int j = 0; j < n; j++) {
            fread(&A[i][j], sizeof(double), 1, f);
        }
    }

    fclose(f);

    Matrix_t eq = {n, A, NULL};

    for(int i = 0; i < 13; i++) {
        double real = IIM(eq, sol, lambda, 1 / pow(10, i), r, 0);
        fprintf(r, "%.16lf;%.16lf\n", 1 / pow(10, i), real);
    }

    fclose(r);
}

int main(void) {
    Sol("iim_cond1_SRC.bin", "iim_cond1_RES.csv", 1e-14);
    Sol("iim_cond2_SRC.bin", "iim_cond2_RES.csv", 1e-14);
    Sol("iim_cond3_SRC.bin", "iim_cond3_RES.csv", 1e-14);
    Sol("iim_cond4_SRC.bin", "iim_cond4_RES.csv", 1e-14);
    Sol("iim_cond5_SRC.bin", "iim_cond5_RES.csv", 1e-14);
    Sol("iim_cond6_SRC.bin", "iim_cond6_RES.csv", 1e-14);

    Conv("iim_cond1_SRC.bin", "iim_cond1_CONV_RES.csv");
    Conv("iim_cond4_SRC.bin", "iim_cond4_CONV_RES.csv");

    return 0;
}
