#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

double** allocate2D(int n, int m) 
{
    double** array = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) 
    {
        array[i] = (double*)malloc(m * sizeof(double));
    }
    return array;
}

double*** allocate3D(int x, int y, int z) 
{
    double*** array = (double***)malloc(x * sizeof(double**));
    for (int i = 0; i < x; i++) 
    {
        array[i] = (double**)malloc(y * sizeof(double*));
        for (int j = 0; j < y; j++) 
        {
            array[i][j] = (double*)malloc(z * sizeof(double));
        }
    }
    return array;
}

complex double** allocate2Dcomplex(int n, int m) 
{
    complex double** array = (complex double**)malloc(n * sizeof(complex double*));
    for (int i = 0; i < n; i++) 
    {
        array[i] = (complex double*)malloc(m * sizeof(complex double));
    }
    return array;
}

complex double*** allocate3Dcomplex(int x, int y, int z) 
{
    complex double*** array = (complex double***)malloc(x * sizeof(complex double**));
    for (int i = 0; i < x; i++) 
    {
        array[i] = (complex double**)malloc(y * sizeof(complex double*));
        for (int j = 0; j < y; j++) 
        {
            array[i][j] = (complex double*)malloc(z * sizeof(complex double));
        }
    }
    return array;
}

void dft(double *input, complex double *output, int N) 
{
    for (int i = 0; i < N; ++i) 
    {
        output[i] = 0.0 + 0.0 * I;
        for (int n = 0; n < N; ++n) 
        {
            double angle = -2.0 * PI * i * n / N;
            output[i] += input[n] * (cos(angle) + I * sin(angle));
        }
    }
}


void dft2d(double **input, complex double **output, int M, int N) {
    for (int u = 0; u < M; ++u) 
    {
        for (int v = 0; v < N; ++v) 
        {
            output[u][v] = 0.0 + 0.0 * I;
            for (int i = 0; i < M; ++i) 
            {
                for (int j = 0; j < N; ++j) 
                {
                    double angle = -2.0 * PI * ((u * i / (double)M) + (v * j / (double)N));
                    output[u][v] += input[i][j] * (cos(angle) + I * sin(angle));
                }
            }
        }
    }
}

int main() 
{
    int N = 30;
    double** input = allocate2D(N, N);
    //Open the file and write the header
    FILE * fp = fopen("input.txt","wb");

    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++)
        {
            input[i][j] = 0;
        }
    }

    input[0][5] = 1;
    input[0][25] = 1;

    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(fp, "%f ", input[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

     // Массив для результата
    complex double** output = allocate2Dcomplex(N, N);

    // Выполнение ДПФ
    dft2d(input, output, N, N);

    fp = fopen("output.txt","wb");
    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(fp, "%f ", creal(output[i][j]));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 0;
}
