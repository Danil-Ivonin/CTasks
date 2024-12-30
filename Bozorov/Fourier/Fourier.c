#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>

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

void clear2d(double** array, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(array[i]);
    }
    free(array);
}

void clear3d(double*** array, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
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


void dft2d(double **input, double *output, int start, int end, int N) 
{ 
    int idx = 0;
    for (int u = start; u <= end; ++u) 
    {
        for (int v = 0; v < N; ++v) 
        {
            complex double tmp = 0.0 + 0.0 * I;
            for (int i = 0; i < N; ++i) 
            {
                for (int j = 0; j < N; ++j) 
                {
                    double angle = -2.0 * PI * ((u * i / (double)N) + (v * j / (double)N));
                    tmp += input[i][j] * (cos(angle) + I * sin(angle));
                }
            }
            output[idx] = creal(tmp);
            idx++;
        }
    }
}

void dft3d(double ***input, double *output, int start, int end, int N) 
{ 
    int idx = 0;
    for (int u = start; u <= end; ++u) 
    {
        for (int v = 0; v < N; ++v) 
        {
            for (int z = 0; z < N; ++z) 
            {
                complex double tmp = 0.0 + 0.0 * I;
                for (int i = 0; i < N; ++i) 
                {
                    for (int j = 0; j < N; ++j) 
                    {
                        for (int k = 0; k < N; ++k) 
                        {
                            double angle = -2.0 * PI * ((u * i / (double)N) + (v * j / (double)N) + (z * k / (double)N));
                            tmp += input[i][j][k] * (cos(angle) + I * sin(angle));
                        }
                    }
                }
                output[idx] = creal(tmp);
                idx++;
            }
        }
    }
}

int main(int argc, char* argv[]) 
{
    int rank, size;
    double tstart, tend;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = 30;
    double** input2d = allocate2D(N, N);
    double*** input3d = allocate3D(N, N, N);

    //fill array
    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++)
        {
            input2d[i][j] = 0;
            for (int k = 0; k < N; k++)
            {
                input3d[i][j][k] = 0;
            }
        }
    }

    input2d[0][5] = 1;
    input2d[0][25] = 1;

    input3d[0][0][5] = 1;
    input3d[0][0][25] = 1;

    if (rank == 0)
    {
        FILE *fp = fopen("input2d.txt","wb");
        FILE *fp2 = fopen("input3d.txt","wb");
        //write data to file
        for (int i = 0; i < N; i++) 
        {
            for (int j = 0; j < N; j++)
            {
                fprintf(fp, "%f ", input2d[i][j]);
                for (int k = 0; k < N; k++)
                {
                    fprintf(fp2, "%f ", input3d[i][j][k]);
                }
            }
            fprintf(fp, "\n");
            fprintf(fp2, "\n");
        }

        fclose(fp);
        fclose(fp2);
    }

    int chunk_size = ceil((double)N / size);
    int start = rank * chunk_size;
    int end = (rank+1) * chunk_size - 1;
    if (end > N)
        end = N - 1;

    double* res2d = (double*)malloc(N * chunk_size * sizeof(double));
    double* res3d = (double*)malloc(N * N * chunk_size * sizeof(double));
    MPI_Barrier(MPI_COMM_WORLD);

    double start2d = MPI_Wtime();
    dft2d(input2d, res2d, start, end, N);
    double end2d = MPI_Wtime();

    double start3d = MPI_Wtime();
    dft3d(input3d, res3d, start, end, N);
    double end3d = MPI_Wtime();

    double* output2d;
    double* output3d;
    if (rank == 0)
    {
        output2d = (double*)malloc(N * chunk_size * size * sizeof(double));
        output3d = (double*)malloc(N * N * chunk_size * size * sizeof(double));
    }

    MPI_Gather(res2d, N * chunk_size, MPI_DOUBLE, output2d, N * chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(res3d, N * N * chunk_size, MPI_DOUBLE, output3d, N * N * chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (rank == 0)
    {
        FILE *fp = fopen("output2d.txt","wb");
        FILE *fp2 = fopen("output3d.txt","wb");
        int idx2d = 0;
        int idx3d = 0;
        for (int i = 0; i < N; i++) 
        {
            for (int j = 0; j < N; j++)
            {
                fprintf(fp, "%f ", output2d[idx2d]);
                idx2d++;
                for (int k = 0; k < N; k++)
                {
                    fprintf(fp2, "%f ", output3d[idx3d]);
                    idx3d++;
                }
            }
            fprintf(fp, "\n");
            fprintf(fp2, "\n");
        }
        fclose(fp);
        fclose(fp2);

        printf("max idx2d = %d\n", idx2d);
        printf("max idx3d = %d\n", idx3d);

        double time2d = end2d - start2d;
        printf("Time taken parallel 2d = %f\n", time2d);

        double time3d = end3d - start3d;
        printf("Time taken parallel 3d = %f\n", time3d);
    }

    clear2d(input2d, N);
    clear3d(input3d, N);
    free(output2d);
    free(output3d);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
