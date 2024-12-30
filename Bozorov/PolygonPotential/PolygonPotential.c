#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define N 10000000  // Number of steps for numerical integration
#define PI 3.141592653589793

// Function to calculate the contribution of potential from one side
double potential_of_side(double x, double y, double z, 
                         double x1, double x2, double y1, double y2, double z1, 
                         double lambda) {
    double potential = 0.0;
    double dx = (x2 - x1) / N;
    double dy = (y2 - y1) / N;

    for (int i = 0; i < N; i++) {
        double x_line = x1 + i * dx;
        double y_line = y1 + i * dy;
        double z_line = z1;
        double r = sqrt(pow(x - x_line, 2) + pow(y - y_line, 2) + pow(z - z_line, 2));
        if (r > 1e-9) { // Avoid division by 0
            potential += lambda / r * sqrt(pow(dx, 2) + pow(dy, 2));
        }
    }
    return potential;
}

// Function to calculate the potential from a polygon frame
double potential_polygon(double x, double y, double z, double radius, double lambda, int sides, int rank, int size) {
    double total_potential = 0.0;

    // Calculate vertices of the polygon
    double angle_step = 2 * PI / sides;
    double vertices[sides][2];
    for (int i = 0; i < sides; i++) {
        vertices[i][0] = radius * cos(i * angle_step);
        vertices[i][1] = radius * sin(i * angle_step);
    }

    // Distribute side computations among processes
    for (int i = rank; i < sides; i += size) {
        int next = (i + 1) % sides;
        double x1 = vertices[i][0];
        double y1 = vertices[i][1];
        double x2 = vertices[next][0];
        double y2 = vertices[next][1];
        total_potential += potential_of_side(x, y, z, x1, x2, y1, y2, 0, lambda);
    }

    return total_potential;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <number_of_sides> <radius>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int sides = atoi(argv[1]);
    double radius = atof(argv[2]);
    double lambda = 100; // Linear charge density (C/m)
    double x = 0.2, y = 0.2, z = 1.0; // Observation point coordinates (m)

    if (sides < 3) {
        if (rank == 0) {
            fprintf(stderr, "Error: Number of sides must be at least 3.\n");
        }
        MPI_Finalize();
        return 1;
    }

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes before timing
    double start_time = MPI_Wtime();

    // Each process computes its contribution to the potential
    double local_potential = potential_polygon(x, y, z, radius, lambda, sides, rank, size);

    // Sum up the results from all processes
    double total_potential = 0.0;
    MPI_Reduce(&local_potential, &total_potential, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes after computations
    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Электростатический потенциал: %.3e V\n", total_potential);
        printf("Время выполнения: %.6f секунд\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
