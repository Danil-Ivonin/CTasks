#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void bubbleSort(int* arr, int n) {
    for (int i = 0; i < n - 1; i++) {
        int swapped = 0;
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                int temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
                swapped = 1;
            }
        }
        if (!swapped) break;
    }
}

void printArray(int* arr, int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

void merge(int* arr, int start, int mid, int end) {
    int n1 = mid - start + 1;
    int n2 = end - mid;

    int* left = (int*)malloc(n1 * sizeof(int));
    int* right = (int*)malloc(n2 * sizeof(int));

    for (int i = 0; i < n1; i++) left[i] = arr[start + i];
    for (int i = 0; i < n2; i++) right[i] = arr[mid + 1 + i];

    int i = 0, j = 0, k = start;
    while (i < n1 && j < n2) {
        if (left[i] <= right[j]) {
            arr[k++] = left[i++];
        } else {
            arr[k++] = right[j++];
        }
    }

    while (i < n1) arr[k++] = left[i++];
    while (j < n2) arr[k++] = right[j++];

    free(left);
    free(right);
}

void mergeSortedChunks(int* arr, int* displs, int* send_counts, int size) {
    for (int i = 1; i < size; i++) {
        merge(arr, 0, displs[i] - 1, displs[i] + send_counts[i] - 1);
    }
}

void writeArrayToFile(const char* filename, int* arr, int n) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++) {
        fprintf(file, "%d ", arr[i]);
    }
    fclose(file);
}

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int arr_size = atoi(argv[1]); 
    int* arr = NULL;
    int* send_counts = NULL;
    int* displs = NULL;

    if (rank == 0) {
        arr = (int*)malloc(arr_size * sizeof(int));
        for (int i = 0; i < arr_size; i++) {
            arr[i] = rand() % 100; 
        }

        writeArrayToFile("original.txt", arr, arr_size);

        send_counts = (int*)malloc(size * sizeof(int));
        displs = (int*)malloc(size * sizeof(int));

        int base_size = arr_size / size;
        int extra = arr_size % size;

        for (int i = 0; i < size; i++) {
            send_counts[i] = (i < extra) ? (base_size + 1) : base_size;
            displs[i] = (i == 0) ? 0 : (displs[i - 1] + send_counts[i - 1]);
        }
    }

    double start_time = MPI_Wtime();

    int local_size = (rank < (arr_size % size)) ? ((arr_size / size) + 1) : (arr_size / size);

    int* local_array = (int*)malloc(local_size * sizeof(int));

    MPI_Scatterv(arr, send_counts, displs, MPI_INT, local_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    bubbleSort(local_array, local_size);

    int* gathered_arr = NULL;
    if (rank == 0) {
        gathered_arr = (int*)malloc(arr_size * sizeof(int));
    }

    MPI_Gatherv(local_array, local_size, MPI_INT, gathered_arr, send_counts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        mergeSortedChunks(gathered_arr, displs, send_counts, size);
        

        writeArrayToFile("sorted.txt", gathered_arr, arr_size);

        free(arr);
        free(send_counts);
        free(displs);
        free(gathered_arr);


    }
    free(local_array);
    double end_time = MPI_Wtime();
    MPI_Finalize();
    if (rank == 0) {
        printf("\nTime measured: %f\n", end_time - start_time);
    }
    return 0;
}
