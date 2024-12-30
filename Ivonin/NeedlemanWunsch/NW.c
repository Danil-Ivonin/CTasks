#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MATCH 1
#define MISMATCH -1
#define GAP -1

// Функция для вычисления максимума из трех чисел
int max(int a, int b, int c)
{
    return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

void read_sequences(const char* file_name, char** seq1, char** seq2)
{
    // Открытие файла
    FILE *file = fopen(file_name, "r");
    if (file == NULL)
    {
        perror("Ошибка открытия файла");
        exit(EXIT_FAILURE);
    }

    // Считывание длины последовательностей
    int seq_length;
    if (fscanf(file, "%d", &seq_length) != 1)
    {
        perror("Ошибка чтения длины последовательностей");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Выделение памяти для seq1 и seq2
    *seq1 = (char*)malloc((seq_length + 1) * sizeof(char)); // +1 для '\0'
    *seq2 = (char*)malloc((seq_length + 1) * sizeof(char));
    if (*seq1 == NULL || *seq2 == NULL)
    {
        perror("Ошибка выделения памяти");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Считывание последовательностей
    if (fscanf(file, "%s", *seq1) != 1 || fscanf(file, "%s", *seq2) != 1)
    {
        perror("Ошибка чтения последовательностей");
        free(*seq1);
        free(*seq2);
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fclose(file);
}

int main(int argc, char** argv) {
    int rank, size;

    // Инициализация MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char* seq1 = NULL;
    char* seq2 = NULL;
    const char* file_name = "sequences.txt";

    // Считываем последовательности из файла
    read_sequences(file_name, &seq1, &seq2);

    int n = strlen(seq1);
    int m = strlen(seq2);

    // Создаем матрицу для локальных вычислений
    int** matrix = (int**)malloc((n + 1) * sizeof(int*));
    for (int i = 0; i <= n; i++)
    {
        matrix[i] = (int*)malloc((m + 1) * sizeof(int));
    }

    // Инициализация первой строки и столбца
    for (int i = 0; i <= n; i++)
    {
        matrix[i][0] = i * GAP;
    }
    for (int j = 0; j <= m; j++)
    {
        matrix[0][j] = j * GAP;
    }

    int* recvcounts = (int*)malloc(size * sizeof(int));
    int* displs = (int*)malloc(size * sizeof(int));

    // Распределение диагоналей между процессами
    double tstart = MPI_Wtime();
    int num_diagonals = n + m;
    for (int diag = 1; diag < num_diagonals; diag++)
    {
        int start_i = (diag <= n) ? diag : n;
        int start_j = (diag <= n) ? 1 : diag - n + 1;
        int num_elements = (diag <= n) ? diag : num_diagonals - diag;

        int chunk_size = (num_elements + size - 1) / size; // деление диагонали
        int local_start = rank * chunk_size;
        int local_end = local_start + chunk_size - 1;
        if (local_end >= num_elements) local_end = num_elements - 1;

        int calc_count = local_end >= local_start ? local_end - local_start + 1 : 0;

        int* recv_buf = (int*)malloc(num_elements * sizeof(int));
        int* send_buf = (calc_count > 0) ? (int*)malloc(calc_count * sizeof(int)) : NULL;

        for (int i = 0; i < size; i++)
        {
            int p_start = i * chunk_size;
            int p_end = p_start + chunk_size - 1;
            if (p_end >= num_elements)
            {
                p_end = num_elements - 1;
            }
            int count = (p_end >= p_start) ? p_end - p_start + 1 : 0;
            recvcounts[i] = count;
            displs[i] = (i == 0) ? 0 : displs[i - 1] + recvcounts[i - 1];
        }

        // Вычисляем локальную часть диагонали
        for (int k = local_start; k <= local_end; k++) {
            int i = start_i - k;
            int j = start_j + k;

            // Проверяем границы массива
            if (i < 1 || j < 1 || i > n || j > m) continue;

            int match_score = (seq1[i - 1] == seq2[j - 1]) ? MATCH : MISMATCH;
            send_buf[k - local_start] = max(
                matrix[i - 1][j - 1] + match_score,
                matrix[i - 1][j] + GAP,
                matrix[i][j - 1] + GAP
            );
        }

        // Синхронизация между процессами
        MPI_Allgatherv(
            send_buf, calc_count, MPI_INT,
            recv_buf, recvcounts, displs, MPI_INT,
            MPI_COMM_WORLD
        );

        // Заполняем локальную часть матрицы
        for (int k = 0; k < num_elements; k++) {
            int i = start_i - k;
            int j = start_j + k;
            if (i < 1 || j < 1 || i > n || j > m) continue;
            matrix[i][j] = recv_buf[k];
        }

        // Освобождаем временные буферы
        free(recv_buf);
        if (send_buf) free(send_buf);
    }
    double tend = MPI_Wtime();


    if (rank == 0)
    {
        printf("Итоговая оценка = %d\n", matrix[n][m]);
        printf("Время = %f\nДля последовательности длиной %d\n", tend - tstart, n);
    }

    // Очистка памяти
    for (int i = 0; i <= n; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(recvcounts);
    free(displs);

    // Завершение MPI
    MPI_Finalize();
    return 0;
}
