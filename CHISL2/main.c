#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define GAUSS_OK 1
#define GAUSS_ERROR 0

int gauss(double** mat, double* b, int n, double** res, int* actionsCnt) {
    *res = malloc(sizeof(double*) * n);
    *actionsCnt = 0;

    //forward
    for(int i = 0; i < n; i++) {
        /*
        if(mat[i][i] == 0) {
            char f = 0;

            for(int j = i; j < n; j++) {
                if(fabs(mat[j][i]) > fabs(mat[i][i])) {
                    double* tmp = mat[i];
                    mat[i] = mat[j];
                    mat[j] = tmp;
                    f = 1;
                }
            }

            if(!f)
                return GAUSS_ERROR;
        }
        */

        for(int j = i + 1; j < n; j++) {
            double m = mat[j][i] / mat[i][i];
            (*actionsCnt)++;
            //printf("m: %lf / %lf = %lf\n", mat[j][i], mat[i][i], m);

            for(int k = i; k < n; k++) {
                mat[j][k] -= m * mat[i][k];
                (*actionsCnt)++;
                //printf("%lf ", mat[j][k]);
            }
            b[j] -= b[i] * m;
            (*actionsCnt)++;
            //printf("%lf\n", b[j]);
        }
    }

    //backward
    for(int i = n - 1; i >= 0; i--) {
        (*res)[i] = b[i];
        (*actionsCnt)++;
        for(int j = i + 1; j < n; j++) {
            (*res)[i] -= mat[i][j] * (*res)[j];
            (*actionsCnt)++;
        }
        (*res)[i] = (*res)[i] / mat[i][i];
        (*actionsCnt)++;
    }

    return GAUSS_OK;
}

void procces_err_to_cond(char* in, char* out) {
    FILE* src = fopen(in, "r");
    FILE* res = fopen(out, "w");
    int cnt;
    int tmp;
    fscanf(src,"%i", &cnt);

    for(int i = 0; i < cnt; i++) {
        int n;
        double cond;
        double* X;
        fscanf(src,"%i", &n);
        fscanf(src,"%lf", &cond);

        double* B = malloc(sizeof(double) * n);
        double** mat = malloc(sizeof(double*) * n);

        for(int j = 0; j < n; j++) {
            mat[j] = malloc(sizeof(double) * n);
            for(int k = 0; k < n; k++) {
                fscanf(src,"%lf", &mat[j][k]);
            }
            fscanf(src,"%lf", &B[j]);
        }

        int code = gauss(mat, B, n, &X, &tmp);

        if(code == GAUSS_OK) {
            fprintf(res, "%lf;", cond);
            for(int j = 0; j < n; j++)
                fprintf(res, "%.16lf;", X[j]);
            fprintf(res, "\n");
        }

        free(X);
        for(int j = 0; j < n; j++)
            free(mat[j]);
        free(mat);
        free(B);
    }
    fclose(src);
    fclose(res);
}

void process_time_to_N(char* in, char* out) {
    FILE* src = fopen(in, "r");
    FILE* res = fopen(out, "w");

    int cnt;
    int tmp;
    fscanf(src,"%i", &cnt);

    for(int i = 0; i < cnt; i++) {
        int n;
        double cond;
        double* X;
        fscanf(src,"%i", &n);
        fscanf(src,"%lf", &cond);

        double* B = malloc(sizeof(double) * n);
        double** mat = malloc(sizeof(double*) * n);

        for(int j = 0; j < n; j++) {
            mat[j] = malloc(sizeof(double) * n);
            for(int k = 0; k < n; k++) {
                fscanf(src,"%lf", &mat[j][k]);
            }
            fscanf(src,"%lf", &B[j]);
        }

        int code = gauss(mat, B, n, &X, &tmp);

        if(code == GAUSS_OK) {
            fprintf(res, "%i;%i\n", n, tmp);
        }

        free(X);
        for(int j = 0; j < n; j++)
            free(mat[j]);
        free(mat);
        free(B);
    }
    fclose(src);
    fclose(res);
}

int main(void) {
    //procces_err_to_cond("err_to_cond_SRC.txt", "err_to_cond_RES.csv");
    //procces_err_to_cond("err_to_cond_lc_SRC.txt", "err_to_cond_lc_RES.csv");
    procces_err_to_cond("low_cond_SRC.txt", "low_cond_RES.csv");
    //procces_err_to_cond("high_cond_SRC.txt", "high_cond_RES.csv");
    //process_time_to_N("time_to_N_SRC.txt", "time_to_N_SRC.csv");
    return 0;
}
