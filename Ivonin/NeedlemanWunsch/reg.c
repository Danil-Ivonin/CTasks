#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Функция для умножения двух матриц
void matrix_multiply(double* A, int rowsA, int colsA, 
                     double* B, int rowsB, int colsB, 
                     double* result) 
{
    if (colsA != rowsB)
    {
        fprintf(stderr, "Matrix dimensions do not match for multiplication\n");
        exit(1);
    }

    for (int i = 0; i < rowsA; i++)
    {
        for (int j = 0; j < colsB; j++)
        {
            result[i * colsB + j] = 0;
            for (int k = 0; k < colsA; k++)
            {
                result[i * colsB + j] += A[i * colsA + k] * B[k * colsB + j];
            }
        }
    }
}

// Функция для вычисления обратной матрицы методом Гаусса-Жордана
int inverse_matrix(double* A, int n)
{
    double* augmented = malloc(n * 2 * n * sizeof(double)); // Расширенная матрица (оригинал + единичная)

    // Создание расширенной матрицы
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            augmented[i * 2 * n + j] = A[i * n + j]; // Исходная матрица
        }
        for (int j = n; j < 2 * n; j++)
        {
            augmented[i * 2 * n + j] = (i == j - n) ? 1 : 0; // Единичная матрица
        }
    }

    // Прямой ход (приведение к верхнетреугольному виду)
    for (int i = 0; i < n; i++)
    {
        // Проверка на нулевой элемент на диагонали
        if (fabs(augmented[i * 2 * n + i]) < 1e-9)
        {
            fprintf(stderr, "Matrix is singular and cannot be inverted\n");
            free(augmented);
            return -1;
        }

        // Нормализация текущей строки
        double diag_element = augmented[i * 2 * n + i];
        for (int j = 0; j < 2 * n; j++)
        {
            augmented[i * 2 * n + j] /= diag_element;
        }

        // Обнуление элементов ниже текущей строки
        for (int k = i + 1; k < n; k++)
        {
            double factor = augmented[k * 2 * n + i];
            for (int j = 0; j < 2 * n; j++)
            {
                augmented[k * 2 * n + j] -= factor * augmented[i * 2 * n + j];
            }
        }
    }

    // Обратный ход (приведение к диагональному виду)
    for (int i = n - 1; i >= 0; i--)
    {
        for (int k = i - 1; k >= 0; k--)
        {
            double factor = augmented[k * 2 * n + i];
            for (int j = 0; j < 2 * n; j++) {
                augmented[k * 2 * n + j] -= factor * augmented[i * 2 * n + j];
            }
        }
    }

    // Копирование обратной матрицы из расширенной в исходную
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = augmented[i * 2 * n + j + n];
        }
    }

    free(augmented);
    return 0; // Успешное выполнение
}

// Ridge Regression
void ridge_regression(double* X, double* y, int rows, int cols, double alpha, double* beta)
{
    double* XT = malloc(cols * rows * sizeof(double));  // Транспонированная матрица X
    double* XTX = malloc(cols * cols * sizeof(double)); // XT * X
    double* I = malloc(cols * cols * sizeof(double));   // Единичная матрица
    double* XTy = malloc(cols * sizeof(double));        // XT * y
    double* XTX_alphaI = malloc(cols * cols * sizeof(double)); // XTX + alpha * I

    // Транспонирование матрицы X
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            XT[j * rows + i] = X[i * cols + j];
        }
    }

    // Вычисление XT * X
    matrix_multiply(XT, cols, rows, X, rows, cols, XTX);

    // Создание единичной матрицы и добавление alpha * I
    for (int i = 0; i < cols * cols; i++)
    {
        I[i] = 0;
    }
    for (int i = 0; i < cols; i++)
    {
        I[i * cols + i] = alpha;
    }
    for (int i = 0; i < cols * cols; i++)
    {
        XTX_alphaI[i] = XTX[i] + I[i];
    }

    // Вычисление XT * y
    matrix_multiply(XT, cols, rows, y, rows, 1, XTy);

    // Решение системы (XTX + alpha * I) * beta = XT * y
    // Здесь предполагается использование функции для обратной матрицы
    inverse_matrix(XTX_alphaI, cols * cols); // Реализация inverse_matrix зависит от вашей библиотеки
    matrix_multiply(XTX_alphaI, cols, cols, XTy, cols, 1, beta);

    // Освобождение памяти
    free(XT);
    free(XTX);
    free(I);
    free(XTy);
    free(XTX_alphaI);
}

// LASSO Regression с использованием градиентного спуска
void lasso_regression(double* X, double* y, int rows, int cols, double alpha, int iterations, double lr, double* beta)
{
    for (int i = 0; i < cols; i++)
    {
        beta[i] = 0; // Инициализация коэффициентов
    }

    for (int iter = 0; iter < iterations; iter++)
    {
        double* gradient = malloc(cols * sizeof(double));

        // Вычисление градиента
        for (int j = 0; j < cols; j++)
        {
            gradient[j] = 0;
            for (int i = 0; i < rows; i++)
            {
                double error = y[i];
                for (int k = 0; k < cols; k++)
                {
                    error -= X[i * cols + k] * beta[k];
                }
                gradient[j] += -2 * X[i * cols + j] * error;
            }
            gradient[j] += alpha * (beta[j] > 0 ? 1 : (beta[j] < 0 ? -1 : 0)); // sign(beta[j])
        }

        // Обновление коэффициентов
        for (int j = 0; j < cols; j++) {
            beta[j] -= lr * gradient[j];
        }

        free(gradient);
    }
}



void read_X_train(const char *filename, double **X_train, int *rows, int *cols)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Ошибка при открытии файла x_train.txt");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%d", rows);
    fscanf(file, "%d", cols);

    // Выделить память для двумерного массива
    *X_train = (double *)malloc(*rows * *cols * sizeof(double));

    // Считать данные
    for (int i = 0; i < *rows * *cols; i++)
    {
        fscanf(file, "%lf", &((*X_train)[i]));
    }

    fclose(file);
}

void read_y_train(const char *filename, double **y_train, int *rows)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Ошибка при открытии файла y_train.txt");
        exit(EXIT_FAILURE);
    }

    // Считать количество строк
    fscanf(file, "%d", rows);

    // Выделить память для одномерного массива
    *y_train = (double *)malloc(*rows * sizeof(double));

    // Считать данные
    for (int i = 0; i < *rows; i++)
    {
        fscanf(file, "%lf", &((*y_train)[i]));
    }

    fclose(file);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double *X_train;
    double *y_train;
    int rows, cols;

    // Считать X_train
    read_X_train("x_train.txt", &X_train, &rows, &cols);

    // Считать y_train
    read_y_train("y_train.txt", &y_train, &rows);

    double* ridge_res = (double*)malloc(cols * sizeof(double));
    double* lasso_res = (double*)malloc(cols * sizeof(double));

    ridge_regression(X_train, y_train, rows, cols, 1.0, ridge_res);
    lasso_regression(X_train, y_train, rows, cols, 0.1, 1000, 0.0001, lasso_res);

    printf("\n");

    for (int i = 0; i < cols; i++)
    {
        printf("%f ", ridge_res[i]);
    }
    printf("\n");
    for (int i = 0; i < cols; i++)
    {
        printf("%f ", lasso_res[i]);
    }
    printf("\n");

    free(X_train);
    free(y_train);
    MPI_Finalize();

    return 0;
}
