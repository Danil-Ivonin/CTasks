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

void bubbleSort(int arr[], int low, int high) 
{
    for (int i = low; i < high - 1; i++) 
    {
        // Last i elements are already in place, so the loop
        // will only num n - i - 1 times
        for (int j = low; j < high - i - 1; j++) 
        {
            if (arr[j] > arr[j + 1])
            {
                swap(&arr[j], &arr[j + 1]);
            }
        }
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

void quickSort(int arr[], int low, int high) {
    if (low < high) {

        // call partition function to find Partition Index
        int pi = partition(arr, low, high);

        // Recursively call quickSort() for left and right
        // half based on Partition Index
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

int main() 
{
    int rank, size, chunk_size;
    double tstart_bubble, tend_bubble, tstart_quick, tend_quick;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int arr_size = 30000;
    int* arr1 = (int*)malloc(size * sizeof(int));
    int* arr2 = (int*)malloc(size * sizeof(int));

    if (rank == 0)
    {
        for (int i = 0; i < size; i++) 
        {
            arr1[i] = rand() % size;
            arr2[i] = rand() % size;
            chunk_size = ceil((double)arr_size / size);
        }
    }

    MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&arr1, MPI_INT, arr_size, 0, MPI_COMM_WORLD);
    MPI_Bcast(&arr2, MPI_INT, arr_size, 0, MPI_COMM_WORLD);

    //-----------------------------------------
    //Start acync part
    MPI_Barrier(MPI_COMM_WORLD);

    int start_i = rank * chunk_size;
    int end_i = (rank+1) * chunk_size;
    if (end_i > arr_size)
        end_i = arr_size - 1;


    int* sub_arr_bubble = (int*)malloc((end_i - start_i) * sizeof(int));
    int* sub_arr_quick = (int*)malloc((end_i - start_i) * sizeof(int));

    int index = 0;
    for (int i = start_i; i < end_i; i++)
    {
        sub_arr_bubble[index] = arr1[i];
        sub_arr_quick[index] = arr2[i];
    }

    tstart_bubble = MPI_Wtime();
    //Bubble sort
    for (int i = 0; i < arr_size - 1; i++) 
    {
    
        for (int j = low; j < high - i - 1; j++) 
        {
            if (arr[j] > arr[j + 1])
            {
                swap(&arr[j], &arr[j + 1]);
            }
        }
    }

    // Calling bubble sort on array arr
    
    quickSort(arr2, 0, arr_size - 1);

    for (int i = 0; i < arr_size; i++)
    {
        printf("%d ", arr1[i]);
    }

    return 0;
}