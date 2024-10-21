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
    /* Parse the command line arguments. */
    if (argc != 7) 
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

    int calc_size = xres * yres;
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

    double x, y; /* Coordinates of the current point in the complex plane. */
    double u, v; /* Coordinates of the iterated point. */
    int i,j; /* Pixel counters */
    int k; /* Iteration counter */   
    for (j = 0; j < yres; j++)
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
                //set blue color
                image[index].blue = k;
                if (k > abs(maxiter / 3))
                {
                    //set red color
                    image[index].red = k - abs(maxiter / 3);
                }
                if(k > abs(maxiter / 2))
                {
                    //set green color
                    image[index].green = k - abs(maxiter / 2);
                }
            };
        }
    }

    /* Open the file and write the header. */
    FILE * fp = fopen("pic.ppm","wb");
    /*write ASCII header to the file*/
    fprintf(fp, "P3\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d %d\n%d\n",
            xmin, xmax, ymin, ymax, maxiter,
            xres, yres, maxiter);

    
    for (int i = 0; i < calc_size; i++)
    {
        fprintf(fp, "%d %d %d\n", image[i].red, image[i].green, image[i].blue);
    }
    fclose(fp);

    return 0;
}