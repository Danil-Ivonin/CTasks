#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

typedef struct
{
    uint8_t red;
    uint8_t green;
    uint8_t blue;
}
pixel_t;

int max(int a, int b) { return a > b ? a : b; }

int min(int a, int b) { return a < b ? a : b; }

int main(int argc, char* argv[])
{
    int rank, size;
    double tstart, tend;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    /* Parse the command line arguments. */
    if (argc != 7 && rank == 0) 
    {
        printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <res>\n", argv[0]);
        printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* The window in the plane. */
    const double xmin = atof(argv[1]);
    const double xmax = atof(argv[2]);
    const double ymin = atof(argv[3]);
    const double ymax = atof(argv[4]);

    /* Maximum number of iterations*/
    const int maxiter = (unsigned short)atoi(argv[5]);

    /* Image size, width is given, height is computed. */
    const int xres = atoi(argv[6]);
    const int yres = (xres*(ymax-ymin))/(xmax-xmin);

    printf("Data:\nxmin = %f, xmax = %f, xres = %d\n", xmin, xmax, xres);
    printf("ymin = %f, ymax = %f, yres = %d\n", ymin, ymax, yres);
    printf("img size = %d\n", (xres * yres));
    printf("max iter = %d\n\n", maxiter);

    int chunk_y = ceil((double)yres / size);
    int start_y = rank * chunk_y;
    int end_y = (rank+1) * chunk_y - 1;
    if (end_y > yres)
        end_y = yres - 1;

    int calc_size = xres * chunk_y;

    printf("chunk_y = %d, calc_size = %d at thred %d\n", chunk_y, calc_size, rank);
    pixel_t* image = (pixel_t *)malloc(calc_size * sizeof(pixel_t));
    for (int i = 0; i < calc_size; i++)
    {
        image[i].red = 0;
        image[i].green = 0;
        image[i].blue = 0;
    }

    /* Precompute pixel width and height. */
    double dx=(xmax-xmin)/xres;
    double dy=(ymax-ymin)/yres;
    printf("Start calc data at %d:\nX : %d -> %d\nY : %d -> %d\nCalc size = %d, Chunk_y = %d\n\ns", 
    rank, 0, xres, start_y, end_y, calc_size, chunk_y);
    double x, y; /* Coordinates of the current point in the complex plane. */
    double u, v; /* Coordinates of the iterated point. */
    int i,j; /* Pixel counters */
    int k; /* Iteration counter */   
    for (j = start_y; j < end_y; j++)
    {
        y = ymax - j * dy;
        for(i = 0; i < xres; i++) 
        {
            double u = 0.0;
            double v = 0.0;
            double u2 = 0;
            double v2 = 0;
            x = xmin + i * dx;
            /* iterate the point */
            for (k = 1; k < maxiter && (u2 + v2 < 4.0); k++) 
            {
                v = 2 * u * v + y;
                u = u2 - v2 + x;
                u2 = u * u;
                v2 = v * v;
            };

            /* compute  pixel color and write it to file */
            if (k < maxiter) 
            {
                int index = (j * xres) + i;
                /* exterior */
                if (k > 0)
                {
                    //set blue color
                    image[index].blue = min(k, 256);
                }
                if (k > 256)
                {
                    //set red color
                    image[index].red = min(k-256, 256);
                }
                if(k > 512)
                {
                    //set green color
                    image[index].green = min(k-512, 256);
                }
            };
        }
    }
    printf("Fill send arr\n");
    int* send_arr = (int*)malloc(calc_size * sizeof(int) * 3);
    printf("Send arr size = %d\n", calc_size * 3);
    for (int i = 0; i < calc_size; i++)
    {
        send_arr[3 * i]     = image[i].red;
        send_arr[3 * i + 1] = image[i].green;
        send_arr[3 * i + 2] = image[i].blue;
    }

    int* recv_arr;
    const int recv_arr_size = xres * yres * 3;
    printf("Recv arr size = %d\n", recv_arr_size);
    if(rank == 0)
    {
        int* recv_arr = (int*)malloc(recv_arr_size * sizeof(int));
    }
    printf("Start mpi gather\n");
    MPI_Gather(send_arr, calc_size * 3, MPI_INT, recv_arr, calc_size * 3, MPI_INT, 0, MPI_COMM_WORLD);
    //recv_arr = send_arr;
    if(rank == 0)
    {
        /* Open the file and write the header. */
        FILE * fp = fopen("pic.ppm","wb");
        /*write ASCII header to the file*/
        fprintf(fp, "P3\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d %d\n%d\n",
                xmin, xmax, ymin, ymax, maxiter,
                xres, yres, 256);

        
        for (int i = 0; i < recv_arr_size; i+=3)
        {
            fprintf(fp, "%d %d %d\n", recv_arr[i], recv_arr[i + 1], recv_arr[i + 2]);
        }
        fclose(fp);
    }
    MPI_Finalize();
    return 0;
}