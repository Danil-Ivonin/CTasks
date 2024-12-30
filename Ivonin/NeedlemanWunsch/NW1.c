#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MATCH 1
#define MISMATCH -1
#define GAP -1

// Функция для вычисления максимума из трех чисел
int max(int a, int b, int c)
{
    return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

// Основная программа
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    // Последовательности для выравнивания
    char seq1[] = "GCT";
    char seq2[] = "GCA";

    int n = strlen(seq1);
    int m = strlen(seq2);

    // Создание и инициализация матрицы оценок
    int** matrix = (int**)malloc((n + 1) * sizeof(int*));
    for (int i = 0; i <= n; i++) {
        matrix[i] = (int*)malloc((m + 1) * sizeof(int));
    }

    // Инициализация первой строки и первого столбца
    for (int i = 0; i <= n; i++) {
        matrix[i][0] = i * GAP;
    }
    for (int j = 0; j <= m; j++) {
        matrix[0][j] = j * GAP;
    }

    // Заполнение матрицы
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int match_score = (seq1[i - 1] == seq2[j - 1]) ? MATCH : MISMATCH;
            matrix[i][j] = max(
                matrix[i - 1][j - 1] + match_score, // Замена или совпадение
                matrix[i - 1][j] + GAP,            // Пропуск в seq2
                matrix[i][j - 1] + GAP             // Пропуск в seq1
            );
        }
    }

    // Обратный проход для восстановления выравнивания
    char aligned_seq1[n + m];
    char aligned_seq2[n + m];
    int index1 = 0, index2 = 0;

    int i = n, j = m;
    while (i > 0 && j > 0) {
        int current_score = matrix[i][j];
        int match_score = (seq1[i - 1] == seq2[j - 1]) ? MATCH : MISMATCH;

        if (current_score == matrix[i - 1][j - 1] + match_score) {
            aligned_seq1[index1++] = seq1[i - 1];
            aligned_seq2[index2++] = seq2[j - 1];
            i--;
            j--;
        } else if (current_score == matrix[i - 1][j] + GAP) {
            aligned_seq1[index1++] = seq1[i - 1];
            aligned_seq2[index2++] = '-';
            i--;
        } else {
            aligned_seq1[index1++] = '-';
            aligned_seq2[index2++] = seq2[j - 1];
            j--;
        }
    }

    // Добавляем оставшиеся символы
    while (i > 0) {
        aligned_seq1[index1++] = seq1[i - 1];
        aligned_seq2[index2++] = '-';
        i--;
    }
    while (j > 0) {
        aligned_seq1[index1++] = '-';
        aligned_seq2[index2++] = seq2[j - 1];
        j--;
    }

    // Разворачиваем строки
    aligned_seq1[index1] = '\0';
    aligned_seq2[index2] = '\0';
    for (int k = 0; k < index1 / 2; k++) {
        char temp = aligned_seq1[k];
        aligned_seq1[k] = aligned_seq1[index1 - k - 1];
        aligned_seq1[index1 - k - 1] = temp;
    }
    for (int k = 0; k < index2 / 2; k++) {
        char temp = aligned_seq2[k];
        aligned_seq2[k] = aligned_seq2[index2 - k - 1];
        aligned_seq2[index2 - k - 1] = temp;
    }

    // Вывод результата
    printf("Выравнивание:\n%s\n%s\n", aligned_seq1, aligned_seq2);

    printf("Матрица выравнивания:\n");
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= m; j++)
        {
            printf("%4d", matrix[i][j]);
        }
        printf("\n");
    }
    // Очистка памяти
    for (int i = 0; i <= n; i++) {
        free(matrix[i]);
    }
    free(matrix);
    MPI_Finalize();

    return 0;
}
