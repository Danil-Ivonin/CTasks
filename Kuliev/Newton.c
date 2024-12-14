#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>

#define WIDTH 1000
#define HEIGHT 1000
#define MAX_ITER 100

typedef struct
{
    int red;
    int green;
    int blue;
}
pixel_t;


// Функция для вычисления корней полинома z^5 - 1 и производной
void newton_method(complex double z, complex double* roots, int num_roots, complex double* result) 
{
    complex double f_z = cpow(z, 5) - 1;  // z^5 - 1
    complex double f_prime_z = 5 * cpow(z, 4);  // Производная 5 * z^4
    
    *result = z - f_z / f_prime_z;
}

// Проверка, к какому из корней приближается точка
int find_root_index(complex double z, complex double* roots, int num_roots) 
{
    double min_distance = 1e10;
    int index = 0;
    
    for (int i = 0; i < num_roots; i++) 
    {
        double distance = cabs(z - roots[i]);
        if (distance < min_distance) 
        {
            min_distance = distance;
            index = i;
        }
    }
    return index;
}

void set_colors(pixel_t* out_colors)
{
    pixel_t red;
    red.red = 255;
    red.green = 0;
    red.blue = 0;
    pixel_t blue;
    blue.red = 0;
    blue.blue = 255;
    blue.green = 0;
    pixel_t green;
    green.red = 0;
    green.blue = 0;
    green.green = 255;
    out_colors[0] = red;
    out_colors[1] = blue;
    out_colors[2] = green;
    out_colors[3] = red;
    out_colors[3].blue = 255;
    out_colors[4] = blue;
    out_colors[4].green = 255;
}

int main(int argc, char* argv[]) 
{
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int chunk_size = ceil((double)HEIGHT / size);
    int start_y = rank * chunk_size;
    int end_y = (rank+1) * chunk_size - 1;
    if (end_y > HEIGHT)
        end_y = HEIGHT - 1;

    int calc_size = WIDTH * (end_y - start_y);


    int num_roots = 5;
    complex double roots[5];
    pixel_t* res = (pixel_t*)malloc(calc_size * sizeof(pixel_t));
    pixel_t* colors = (pixel_t*)malloc(5 * sizeof(pixel_t));
    set_colors(colors);

    for (int k = 0; k < 5; k++) 
    {
        roots[k] = cexp(2 * M_PI * I * k / 5);  // Корни z^5 = 1
    }

    double start_time = MPI_Wtime();
    int i = 0;
    // Основной цикл рендеринга
    for (int y = start_y; y < end_y; y++) 
    {
        for (int x = 0; x < WIDTH; x++) 
        {
            // Преобразование пикселя в комплексную плоскость
            double real = 2.0 * (x - WIDTH / 2) / (WIDTH / 2);
            double imag = 2.0 * (y - HEIGHT / 2) / (HEIGHT / 2);
            complex double z = real + imag * I;

            // Итерации метода Ньютона
            complex double new_z = z;
            for (int iter = 0; iter < MAX_ITER; iter++) 
            {
                newton_method(new_z, roots, num_roots, &new_z);
            }

            // Нахождение ближайшего корня
            int root_index = find_root_index(new_z, roots, num_roots);
            
            res[i] = colors[root_index];
            i++;
        }
    }

    double end_time = MPI_Wtime();

    int* send_arr = (int*)malloc(calc_size * sizeof(int) * 3);

    for (int i = 0; i < calc_size; i++)
    {
        send_arr[3 * i]     = res[i].red;
        send_arr[3 * i + 1] = res[i].green;
        send_arr[3 * i + 2] = res[i].blue;
    }

    int* recv_arr;
    const int recv_arr_size = WIDTH * HEIGHT * 3;
    if(rank == 0)
    {
        recv_arr = (int*)malloc(recv_arr_size * sizeof(int));
    }

    MPI_Gather(send_arr, calc_size * 3, MPI_INT, recv_arr, calc_size * 3, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        //Open the file and write the header
        FILE * fp = fopen("pic.ppm","wb");

        //write header to the file
        fprintf(fp, "P3\n# Newton\n%d %d\n%d\n",
                WIDTH, HEIGHT, 255);

        for (int i = 0; i < recv_arr_size; i += 3) 
        {
                fprintf(fp, "%d %d %d\n", recv_arr[i], recv_arr[i + 1], recv_arr[i + 2]);
        }

        fclose(fp);

        double time_taken_parallel = end_time - start_time;
        printf("Time taken parallel = %f for:\nimage size = %dx%d\ncalc depth = %d\n", time_taken_parallel, WIDTH, HEIGHT, MAX_ITER);
    }
    MPI_Finalize();
    return 0;
}
