#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**
 * Вычитает предыдущую строку матрицы из всех последующих строк, которые присутствуют в текущем процессе.
 * @param recv_id       ID процесса, с которого была получена строка,
 * @param id            ID текущего процесса,
 * @param num_eq        Число уравнений в системе,
 * @param proc_rows     Подматрица коэффицентов левой части в текущем процессе,
 * @param proc_vals     Подвектор значений правой части в текущем процессе,
 * @param curr          Строка, которая будет обработана следующей текущим процессом,
 * @param recvd_row     Вся полученная строка,
 * @param rows_per_proc Количество строк в рамках текущего процесса,
 * @param num_proc      Число процессов,
 * @param var_perm  	Перестановки переменных.
 */
void perform_elimination(int recv_id, int id, int num_eq, double* proc_rows, double* proc_vals, int curr, double* recvd_row, int rows_per_proc, int num_proc, int * var_perm) {
	// получение индекса ведущего элемента (pivot) в полученной строке
	int piv_idx = (int) recvd_row[num_eq];
	int tmp = var_perm[piv_idx];
	var_perm[piv_idx] = var_perm[recv_id];
	var_perm[recv_id] = tmp;
	
	for(int i = 0; i < rows_per_proc; i++) {
		// пропускаем строку с pivot-элементом, так как её не нужно обрабатывать
        if(i * num_proc + id == recv_id) {
			continue;
		}
		
        // получение индекса pivot-элемента
        int piv = (int) recvd_row[num_eq];
		// обмен значениями, чтобы учесть перестановку столбцов
        double tmp = proc_rows[(i * num_eq) + recv_id];
		proc_rows[(i * num_eq) + recv_id] = proc_rows[(i * num_eq) + piv];
		proc_rows[(i * num_eq) + piv] = tmp;
	  
        // вычитание pivot-элемента, умноженного на коэффициент для формирование нулей
		if(i >= curr / num_proc) {
			double piv_val = proc_rows[(i*num_eq) + recv_id];
			for(int j = 0; j < num_eq; j++) {
				proc_rows[(i * num_eq) + j] -= (recvd_row[j] * piv_val);
			}
			proc_vals[i] -= recvd_row[num_eq + 1] * piv_val;
		}
	}
}

/**
 * Определение pivot-элемента в строке
 * @param curr      Индекс строки, в которой нужно найти pivot-элемент,
 * @param proc_rows Число строк в текущем процессе,
 * @param num_proc  Общее число процессов,
 * @param num_eq    Число уравнений в системе,
 *
 * @return Индекс pivot-элемента.
 */
int compute_pivot(int curr, int num_proc, int num_eq, double* proc_rows) {
	// переменная для хранения максимального значения в строке
    double mx = -1;
	// индекс pivot-элемента
    int pivot;
	int row_id = curr / num_proc;

	for(int i = curr; i < num_eq; i++) {
		double val = proc_rows[(row_id * num_eq) + i];
		// модуль значения
        if(val < 0)
			val *= -1;
		// определение нового максимума
        if(val > mx) {
			mx = val;
			pivot = i;
		}
	}
	return pivot;
}

/**
 * Деление строки на pivot-элемент
 * @param id        	ID текущего процесса,
 * @param curr      	Индекс строки, в которой будет происходить деление,
 * @param proc_rows 	Строка в текущем процессе,
 * @param pivot     	Индекс pivot-элемента
 * @param num_proc  	Число процессов
 * @param num_eq    	Число уравнения в системе, The number of equations in the given systemc_vals The chink of the values vector contained within the current process
 * @param rows_per_proc	Число строк в процессе
 * @param proc_vals		Правые части уравнений (представлены в виде одного массива)
 */
void perform_division(int id, int curr, double* proc_rows, int pivot, int num_proc, int num_eq, int rows_per_proc, double * proc_vals) {
	int row_id = curr / num_proc;
    // меняем значения в curr и pivot, чтобы выполнить деление
	double tmp = proc_rows[(row_id * num_eq) + pivot]; 
	proc_rows[(row_id * num_eq) + pivot] = proc_rows[(row_id * num_eq) + curr];
	proc_rows[(row_id * num_eq) + curr] = tmp;
	
    // получаем значение пивота в текущей строке
	double piv_val = proc_rows[(row_id * num_eq) + curr];
	
    // делим все элементы строки на значение pivot-элемента
    for(int i = curr; i < num_eq; i++) {
		proc_rows[(row_id*num_eq) + i] /= piv_val;
	}
    // делим все правые части на pivot-элемент
	proc_vals[row_id] /= piv_val;
}

int main(int argc, char** argv) {
	int num_proc;      // число процессов
	int id;            // индекс процесса
	int num_eq = 0;    // число уравнений
	int rows_per_proc; // число строк на каждый процесс
	double time_taken; // затраченное время

	// Инициализация MPI
	MPI_Init(&argc, &argv);

	// Получение числа процессов
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	
	// Получение текущего индекса процесса
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	

	// Определение потока ввода данных из файла
    FILE* inputfp = NULL;
	if(id == 0) {
		inputfp = fopen(argv[1], "r");
		
		/* read number of equations */
		fscanf(inputfp, "%d", &num_eq);
	}

    // Определение массивов для хранения коэффициентов уравнений
	
	double eq_mat[num_eq * num_eq]; // матрица коэффициентов левых частей (представлена в виде 1D массива)
	double val_mat[num_eq]; 	    // вектор правых частей

	int divs[num_proc];   // число коэффициентов на процесс
	int displs[num_proc]; // начальные индексы коэффициентов каждого процесса

    // Начало параллельных вычислений
	
    time_taken = MPI_Wtime();


	MPI_Bcast(&num_eq, 1, MPI_INT, 0, MPI_COMM_WORLD); // Передача числа уравнений на все процессы
	int rpp = num_eq / num_proc; // Количество строк на процесс по умолчанию 
                                 // (в случае равномерного распредения строк между процессами)
	 
	int np0 = num_proc - (num_eq % num_proc); 	// Количество процессов, которым назначено ровно rpp строк
	int np1 = num_proc - np0; 					// Количество процессов, которым назначено rpp + 1 строк.

	// определение количества строк для текущего процесса
	if(id >= np1)
		rows_per_proc = rpp;
	else
		rows_per_proc = rpp + 1;
	
	// определение массива для перестановок (по умолчанию перестановок нет)
	int var_perm[num_eq];
	for(int i = 0; i < num_eq; i++) {
		var_perm[i] = i;
	}
	
	// чтение данных из файла на мастер-процессе
	if(id == 0) {    
		
		// инициализация массива для подсчета количества строк на каждом процессе
        for(int i = 0; i < num_proc; i++) {
			divs[i] = 0;
		}
		
		// распределение строк между процессами
		for(int i = 0; i < num_eq; i++) {
			
            int assigned_proc = i % num_proc; // определение процесса, которому назначена текущая строка
			
			// Вычисление фактической позиции строки в разбиении
			int rows_before = 0; 
			if(assigned_proc <= np1) {
                rows_before = (rpp + 1) * (assigned_proc);
			} else {
				rows_before = (rpp + 1) * (np1);
				rows_before += (assigned_proc - np1)*(rpp);
			}

			int eff_row = rows_before + (divs[assigned_proc]++);
			
			// Чтение коэффициентов для левой и правой части уравнения
			for(int j = 0; j < num_eq; j++) {
				fscanf(inputfp, "%lf", &eq_mat[eff_row * num_eq + j]);
			}
			fscanf(inputfp, "%lf", &val_mat[eff_row]);
		}

        // определения вектора смещений
		displs[0] = 0;
		for(int i = 1; i < num_proc; i++) {
			divs[i - 1] *= num_eq;
			displs[i] = displs[i - 1] + divs[i - 1];
		}
		divs[num_proc - 1] *= num_eq;
		
        // закрытие файла
		fclose(inputfp);
	}

	// Ожидание завершения чтения данных всеми процессами.
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Разбиение данных на чанки (фрагменты) для каждого процесса
	double proc_rows[num_eq * rows_per_proc]; // данные для текущего процесса (матрица коэффициентов)
    // Распределение строк матрицы коэффициентов между процессами.
	MPI_Scatterv(eq_mat, divs, displs, MPI_DOUBLE, proc_rows, rows_per_proc * num_eq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// Разбиение значений правых частей на чанки для каждого процесса
	double proc_vals[rows_per_proc]; // Вектор значений правых частей для текущего процесса
	for(int i = 0; i < num_proc; i++) {
		divs[i] /= num_eq;
		displs[i] /= num_eq;
	}
	MPI_Scatterv(val_mat, divs, displs, MPI_DOUBLE, proc_vals, rows_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// Определение процесса, который будет получать данные от текущего процесса
	int next_proc = (id + 1) % num_proc;

	// Определение процесса, который будет отправлять данные текущему процессу
	int prev_proc = (id - 1 + num_proc) % num_proc;
	
	int curr = id; 					// текущий процесс обрабатывает строку с индексом id.
	int prev_curr = -1; 			// индекс предыдущей строки
	int piv; 						// индекс пивота для текущей строки
	double recvd_row[num_eq + 2]; 	// массив для хранения строки, полученной от другого процесса + пивот и значение правой части
	MPI_Status st; 					// переменная для хранения статуса передачи сообщения

	/*************** Формирование верхней треугольной матрицы ***************/

	// итерация по всем строкам в процессе
    for(int cnt = 0; cnt < rows_per_proc; cnt++) {
		// начиная с предыдщей строки по текущую
		for(int i = prev_curr + 1; i < curr; i++) {
			
            // принимаем коэффициенты левой части, pivot-элемент и правую часть строки
			MPI_Recv(recvd_row, num_eq + 2, MPI_DOUBLE, prev_proc, i, MPI_COMM_WORLD, &st);
			
            // после получения строки от предыдущего процесса, проверяется, не является ли это последней обработанной строкой.
			if(curr - i < num_proc - 1) {
				MPI_Send(recvd_row, num_eq + 2, MPI_DOUBLE, next_proc, i, MPI_COMM_WORLD);
			}

			// вычитание
			perform_elimination(i, id, num_eq, proc_rows, proc_vals, curr, recvd_row, rows_per_proc, num_proc, var_perm);
		}
		// определение индекса pivot-элемента
        piv = compute_pivot(curr, num_proc, num_eq, proc_rows); 

		// деление на pivot-элемент
        perform_division(id, curr, proc_rows, piv, num_proc, num_eq, rows_per_proc, proc_vals);
		
        // Подготовка данных для отправки следующему процессу:
	    // Заполняем буфер данными для отправки (коэффициенты текущей строки, pivot и правую часть)
		double send_buf[num_eq + 2];
		for(int j = 0; j < num_eq; j++) {
			send_buf[j] = proc_rows[(cnt * num_eq) + j];
		}
		send_buf[num_eq] = piv; // Добавляем индекс pivot-элемента
		send_buf[num_eq + 1] = proc_vals[cnt]; // Добавляем правую часть уравнения для текущей строки
		
		// Отправляем данные следующему процессу, если процессов больше одного
		if(num_proc > 1) {
			MPI_Send(send_buf, num_eq + 2, MPI_DOUBLE, next_proc, curr, MPI_COMM_WORLD);
		}

        // Обновляем индексы текущей и предыдущей строки для дальнейшей обработки
		prev_curr = curr;
		curr += num_proc;
		
        // Вычитание строк с учётом полученной строки и перемещения pivot-элементов
		perform_elimination(prev_curr, id, num_eq, proc_rows, proc_vals, curr, send_buf, rows_per_proc, num_proc, var_perm);
	}

	// Обработка оставшихся строк, которые не будут изменяться, но нужны для обмена с другими процессами
	for(int i = prev_curr + 1; i < num_eq; i++) {
		// Получаем строку, которая будет использоваться для обмена
        MPI_Recv(recvd_row, num_eq + 2, MPI_DOUBLE, prev_proc, i, MPI_COMM_WORLD, &st);
		
		// Если текущая строка не последняя для обработки, отправляем её следующему процессу
        if(curr - i < num_proc - 1) { 
			MPI_Send(recvd_row, num_eq + 2, MPI_DOUBLE, next_proc, i, MPI_COMM_WORLD);
		}
        // Выполняем вычитание полученной строки
		perform_elimination(i, id, num_eq, proc_rows, proc_vals, curr, recvd_row, rows_per_proc, num_proc, var_perm);
	}

    // Ожидаем завершения всех операций обмена данными между процессами
	MPI_Barrier(MPI_COMM_WORLD);

	/*************** Back Substitution Phase ***************/
	
    // Вычисление индекса последней строки, которую будет обрабатывать текущий процесс
	int lst = id + (rows_per_proc - 1) * num_proc;
	int prev_lst = num_eq; // переменная для отслеживания предыдущей последней строки
	double x_val; // значение для хранения результата вычислений
	int count = 1; // количество значений, которое нужно передавать через MPI (в данном случае одно)

	// Основной цикл обратного хода для решения системы
    for(int cnt = rows_per_proc - 1; cnt >= 0; cnt--) {
		// Итерация по всем строкам от последней строки до первой в текущем процессе
        for(int i = prev_lst - 1; i > lst; i--) {
			// Получаем значение переменной x из следующего процесса, которое будет использоваться для обратного хода
            MPI_Recv(&x_val, count, MPI_DOUBLE, next_proc, i + num_eq, MPI_COMM_WORLD, &st);
			
			// Если строка ещё не последняя для обработки, отправляем значение x обратно в предыдущий процесс
            if(i - lst < num_proc - 1) {
				MPI_Send(&x_val, count, MPI_DOUBLE, prev_proc, i + num_eq, MPI_COMM_WORLD);
			}
			
            // Вычитание соответствующего элемента из правой части уравнения (процесс обратной подстановки)
			for(int k = cnt; k >= 0; k--) 
				proc_vals[k] -= proc_rows[(k * num_eq) + i] * x_val;
		}
		// Решение для текущей строки после обратного хода
        double ans = proc_vals[cnt];

		// Обновление правой части уравнения для всех предыдущих строк
        for(int k = cnt - 1; k >= 0; k--)
			proc_vals[k] -= proc_rows[(k * num_eq) + lst] * ans;
		
		// Если процесс не является первым, и строки ещё не обработаны, отправляем результат в предыдущий процесс
        if(lst > 0 && num_proc > 1) {
			// Не отправляем первый pivot, потому что он не нужен для следующих процессов
			MPI_Send(&ans, count, MPI_DOUBLE, prev_proc, lst + num_eq, MPI_COMM_WORLD);
		}

        // Обновление индексов последней строки для следующей итерации
		prev_lst = lst;
		lst -= num_proc;
	}

	// Массив для сбора результатов от всех процессов
    double res[num_eq];

	// Сбор результатов из всех процессов с помощью MPI_Gatherv
    MPI_Gatherv(proc_vals, rows_per_proc, MPI_DOUBLE, res, divs, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// Если текущий процесс - главный (id == 0), продолжаем обработку решения
    if(id == 0) {
		
		// Индекс для заполнения массива решения
        int k = 0;
		double solution[num_eq];
		
		// Переставляем данные в последовательный формат (глобальный порядок переменных)
		for(int i = 0; i < num_proc; i++) {
			for(int j = i; j < num_eq; j += num_proc) {
				solution[j] = res[k++];
			}
		}

        // Массив для хранения конечного решения с учетом перестановки переменных
		double sol[num_eq];
		for(int i = 0; i < num_eq; i++) 
			sol[var_perm[i]] = solution[i];
		
		// Выводим решение в файл
		FILE* outfile = fopen("output.txt", "w");
		fprintf(outfile, "The solution of the system of equations is: \n");
		for(int i = 0; i < num_eq; i++) {
			fprintf(outfile, "x%d = %lf \n", i + 1, sol[i]);
		}
		fclose(outfile);
	}
	  
	// Ожидаем завершения всех операций перед выходом, чтобы MPI_Finalize не был вызван раньше времени
	MPI_Barrier(MPI_COMM_WORLD); 

	// Подсчёт времени выполнения программы
    time_taken = MPI_Wtime() - time_taken;

	if(id == 0) {
		printf("Time taken for execution is %lf\n", time_taken);
	}

	MPI_Finalize();
	return 0;
}