#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

typedef struct
{
    int red;
    int green;
    int blue;
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
    
    int maxiter, xres, yres, chunk_size;
    double xmin, xmax, ymin, ymax, dx, dy;

    //------------------------------------
    // Parse the command line arguments
    if (rank == 0) 
    {
        if (argc != 7)
        {
            printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <res>\n", argv[0]);
            printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024\n", argv[0]);
            exit(EXIT_FAILURE);
        }
        
        //Get user input data
        xmin = atof(argv[1]);
        xmax = atof(argv[2]);
        ymin = atof(argv[3]);
        ymax = atof(argv[4]);

        maxiter = (unsigned short)atoi(argv[5]);

        //Image size, width is given, height is computed
        xres = atoi(argv[6]);
        yres = (xres*(ymax-ymin))/(xmax-xmin);


        // Precompute pixel width and height
        dx = (xmax-xmin)/xres;
        dy = (ymax-ymin)/yres;

        printf("Data:\nxmin = %f, xmax = %f, xres = %d\n", xmin, xmax, xres);
        printf("ymin = %f, ymax = %f, yres = %d\n", ymin, ymax, yres);
        printf("dx = %f, dy = %f\n", dx, dy);
        printf("img size = %d\n", (xres * yres));
        printf("max iter = %d\n\n", maxiter);

        chunk_size = ceil((double)yres / size);
    }

    MPI_Bcast(&maxiter, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xres, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&yres, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ymin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ymax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //-----------------------------------------
    //Start acync part
    MPI_Barrier(MPI_COMM_WORLD);

    printf("Data on thread %d:\n", rank);
    printf("maxiter = %d, xres = %d, yres = %d\n", maxiter, xres, yres);
    printf("xmin = %f, xmax = %f, ymin = %f, ymax = %f\n", xmin, xmax, ymin, ymax);

    int start_y = rank * chunk_size;
    int end_y = (rank+1) * chunk_size - 1;
    if (end_y > yres)
        end_y = yres - 1;

    int calc_size = xres * (end_y - start_y + 1);
    printf("chunk_y = %d, calc_size = %d\n\n", chunk_size, calc_size);

    pixel_t* image = (pixel_t *)malloc(calc_size * sizeof(pixel_t));
    //init black image
    for (int i = 0; i < calc_size; i++)
    {
        image[i].red = 0;
        image[i].green = 0;
        image[i].blue = 0;
    }


    printf("Start calc data at %d:\nX : %d -> %d\nY : %d -> %d\nCalc size = %d, Chunk_y = %d\n\n", 
    rank, 0, xres, start_y, end_y, calc_size, chunk_size);
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
                int bias = j > start_y ? j : (j - start_y);
                int index = ((j - start_y) * xres) + i;
                /* exterior */
                image[index].blue = k;
                image[index].green = k;
                image[index].red = k;
                /*
                if (k > abs(maxiter / 3))
                {
                    //set red color
                    image[index].red = k - abs(maxiter / 3);
                }
                if(k > abs(maxiter / 2))
                {
                    //set green color
                    image[index].green = k - abs(maxiter / 2);
                }*/
            };
        }
    }
    printf("Fill send arr\n");
    int* send_arr = (int*)malloc(calc_size * sizeof(int) * 3);
    printf("Send arr size = %d\n\n", calc_size * 3);
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
        recv_arr = (int*)malloc(recv_arr_size * sizeof(int));
    }
    printf("Start mpi gather\n\n");
    MPI_Gather(send_arr, calc_size * 3, MPI_INT, recv_arr, calc_size * 3, MPI_INT, 0, MPI_COMM_WORLD);
    //recv_arr = send_arr;
    if(rank == 0)
    {
        /* Open the file and write the header. */
        FILE * fp = fopen("pic.ppm","wb");
        /*write ASCII header to the file*/
        fprintf(fp, "P3\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d %d\n%d\n",
            xmin, xmax, ymin, ymax, maxiter,
            xres, yres, maxiter);

        
        for (int i = 0; i < recv_arr_size; i+=3)
        {
            fprintf(fp, "%d %d %d\n", recv_arr[i], recv_arr[i + 1], recv_arr[i + 2]);
        }
        fclose(fp);
    }
    MPI_Finalize();
    return 0;
}