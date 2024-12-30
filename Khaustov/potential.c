#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define N 10000000  // Number of steps for numerical integration

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

// Function to calculate the potential from a square frame
double potential_square(double x, double y, double z, double a, double lambda, int rank, int size) {
    double total_potential = 0.0;

    // Coordinates of the square sides (in the plane z = 0)
    double sides[4][4] = {
        {0, a, 0, 0},  // Bottom side
        {a, a, 0, a},  // Right side
        {a, 0, a, a},  // Top side
        {0, 0, a, 0}   // Left side
    };

    // Distribute side computations among processes
    for (int i = rank; i < 4; i += size) {
        double x1 = sides[i][0];
        double x2 = sides[i][1];
        double y1 = sides[i][2];
        double y2 = sides[i][3];
        total_potential += potential_of_side(x, y, z, x1, x2, y1, y2, 0, lambda);
    }

    return total_potential;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double a = 1.0;      // Side length of the square (m)
    double lambda = 100; // Linear charge density (C/m)
    double x = 0.2, y = 0.2, z = 1.0; // Observation point coordinates (m)

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes before timing
    double start_time = MPI_Wtime();

    // Each process computes its contribution to the potential
    double local_potential = potential_square(x, y, z, a, lambda, rank, size);

    // Sum up the results from all processes
    double total_potential = 0.0;
    MPI_Reduce(&local_potential, &total_potential, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes after computations
    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Electrostatic potential at point (%.2f, %.2f, %.2f): %.3e V\n", x, y, z, total_potential);
        printf("Execution time: %.6f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
