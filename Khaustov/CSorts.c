#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void swap(int* a, int* b) 
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

void bubble_sort(int* arr, int arr_len, int rank, int size) 
{
    int chunk_size = ceil((double)arr_len / size);
    int start = rank * chunk_size;
    int end = (rank + 1) * chunk_size - 1;
    if (end > arr_len)
        end = arr_len - 1;

    for (int i = 0; i < arr_len - 1; i++)
    {
        if (rank == 1)
        {
            printf("\n");
        }
        for (int j = start; j < end; j++) 
        {
            if (rank == 1)
            {
                printf("%d %d\t", arr[j], arr[j+1]);
            }
            if (arr[j] > arr[j + 1])
            {
                swap(&arr[j], &arr[j + 1]);
            }
        }
        if (rank == 1)
        {
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int partition(int arr[], int low, int high) 
{
    // Initialize pivot to be the first element
    int p = arr[low];
    int i = low;
    int j = high;

    while (i < j) 
    {

        // Find the first element greater than
        // the pivot (from starting)
        while (arr[i] <= p && i <= high - 1) {
            i++;
        }

        // Find the first element smaller than
        // the pivot (from last)
        while (arr[j] > p && j >= low + 1) {
            j--;
        }
        if (i < j) {
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[low], &arr[j]);
    return j;
}

void quick_sort(int arr[], int low, int high) 
{
    if (low < high)
    {
        printf("\ninput array: ");
        for (int i = 0; i < 10; i++)
        {
            printf("%d", arr[i]);
        }
        int i = low;
        int j = high;
        int p = low + (high - low) / 2;
        int pivot = arr[p];
        printf("\npivot: %d:%d\n", p, pivot);
        while (i < j)
        {
            while(arr[i] <= pivot && i <= high - 1)
            {
                i++;
            }

            while(arr[j] > pivot && j >= low + 1)
            {
                j--;
            }

            if(i < j)
            {
                printf("swap %d:%d -> %d:%d\n", i, arr[i], j, arr[j]);
                swap(&arr[i], &arr[j]);
                for (int k = 0; k < 10; k++)
                {
                    printf("%d", arr[k]);
                }
                printf("\n");
            }
        }

        printf("\nresult array: ");
        for (int i = 0; i < 10; i++)
        {
            printf("%d", arr[i]);
        }
        printf("\n");
        quick_sort(arr, low, j-1);
        quick_sort(arr, j+1, high);
    }
}


int main(int argc, char* argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int arr_size = 10;
    int* arr = (int*)malloc(arr_size * sizeof(int));
    
    if (rank == 0)
    {
        printf("\n");
        for (int i = 0; i < arr_size; i++)
        {
            arr[i] = rand() % arr_size;
            printf("%d", arr[i]);
        }
        printf("\n");
    }

    MPI_Bcast(arr, arr_size, MPI_INT, 0, MPI_COMM_WORLD);

    bubble_sort(arr, arr_size, rank, size);

    if (rank == 1)
    {
        for (int i = 0; i < arr_size; i++)
        {
            printf("%d", arr[i]);
        }
        printf("\n");
    }
    MPI_Finalize();
    return 0;
}