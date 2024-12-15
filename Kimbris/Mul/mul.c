#include <stdio.h>
#include <stdlib.h>

void readMatrixFromFile(const char *filename, int **matrix, int n) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Unable to open the file");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%d", &matrix[i][j]);
        }
    }

    fclose(file);
}

void multiplyMatrices(int **matrix1, int **matrix2, int **result, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = 0;
            for (int k = 0; k < n; k++) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}

void printMatrix(int **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <size of matrix (n)> <path to matrix file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Get matrix size from command line argument
    int n = atoi(argv[1]);

    // Validate matrix size
    if (n <= 0) {
        fprintf(stderr, "Matrix size must be a positive integer.\n");
        return EXIT_FAILURE;
    }

    // Get file path from command line argument
    const char *filePath = argv[2];

    // Dynamically allocate memory for the matrices
    int **matrix1 = (int **)malloc(n * sizeof(int *));
    int **matrix2 = (int **)malloc(n * sizeof(int *));
    int **result = (int **)malloc(n * sizeof(int *));

    for (int i = 0; i < n; i++) {
        matrix1[i] = (int *)malloc(n * sizeof(int));
        matrix2[i] = (int *)malloc(n * sizeof(int));
        result[i] = (int *)malloc(n * sizeof(int));
    }

    // Read the matrix from the specified file
    readMatrixFromFile(filePath, matrix1, n);

    // Matrix1 and Matrix2 are the same for multiplication
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix2[i][j] = matrix1[i][j]; // Copy matrix1 to matrix2
        }
    }

    // Multiply the matrix with itself
    multiplyMatrices(matrix1, matrix2, result, n);

    // Free dynamically allocated memory
    for (int i = 0; i < n; i++) {
        free(matrix1[i]);
        free(matrix2[i]);
        free(result[i]);
    }
    free(matrix1);
    free(matrix2);
    free(result);

    return 0;
}
