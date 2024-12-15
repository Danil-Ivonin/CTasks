#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void read_matrix_from_file(const char* filename, int* matrix, int size) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < size * size; i++) {
        if (fscanf(file, "%d", &matrix[i]) != 1) {
            perror("Error reading matrix");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    fclose(file);
}

void print_matrix(int* matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%d ", matrix[i * size + j]);
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) {
            printf("Usage: mpirun -np <num_procs> ./mpi_mul <matrix_file> <matrix_size>\n");
        }
        MPI_Finalize();
        return 1;
    }

    const char* filename = argv[1];
    int N = atoi(argv[2]);
    
    if (N % size != 0) {
        if (rank == 0) {
            printf("Matrix size should be divisible by the number of processes\n");
        }
        MPI_Finalize();
        return 1;
    }

    int* A = NULL;
    int* B = NULL;
    int* C = NULL;
    int local_N = N / size;

    if (rank == 0) {
        A = (int*)malloc(N * N * sizeof(int));
        B = (int*)malloc(N * N * sizeof(int));
        C = (int*)malloc(N * N * sizeof(int));
        read_matrix_from_file(filename, A, N);

        // Copy matrix A into B (since we are multiplying A * A)
        for (int i = 0; i < N * N; i++) {
            B[i] = A[i];
        }
    }

    int* local_A = (int*)malloc(local_N * N * sizeof(int));
    int* local_C = (int*)malloc(local_N * N * sizeof(int));

    // Scatter the rows of matrix A to all processes
    MPI_Scatter(A, local_N * N, MPI_INT, local_A, local_N * N, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast matrix B to all processes (all processes need the entire matrix B)
    if (rank != 0) {
        B = (int*)malloc(N * N * sizeof(int));
    }
    MPI_Bcast(B, N * N, MPI_INT, 0, MPI_COMM_WORLD);

    // Perform local matrix multiplication
    for (int i = 0; i < local_N; i++) {
        for (int j = 0; j < N; j++) {
            local_C[i * N + j] = 0;
            for (int k = 0; k < N; k++) {
                local_C[i * N + j] += local_A[i * N + k] * B[k * N + j];
            }
        }
    }

    // Gather the results back to process 0
    MPI_Gather(local_C, local_N * N, MPI_INT, C, local_N * N, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        free(A);
        free(B);
        free(C);
    }

    free(local_A);
    free(local_C);

    MPI_Finalize();
    return 0;
}
