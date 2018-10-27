#include <stdio.h>
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <chrono>
#include <iostream>

using namespace std::chrono;
using namespace std;

// количество строк в исходной квадратной матрице
const int MATRIX_SIZE = 1500;


/// Функция InitMatrix() заполняет переданную в качестве 
/// параметра квадратную матрицу случайными значениями
/// matrix - исходная матрица СЛАУ
void InitMatrix(double** matrix)
{
	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		matrix[i] = new double[MATRIX_SIZE + 1];
	}

	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		for (int j = 0; j <= MATRIX_SIZE; ++j)
		{
			matrix[i][j] = rand() % 2500 + 1;
		}
	}
}

/// Функция SerialGaussMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбец матрицы - значения правых частей уравнений
/// rows - количество строк в исходной матрице
/// result - массив ответов СЛАУ
duration<double> SerialGaussMethod_posled(double **matrix, const int rows, double* result)
{
	int k;
	double koef;

	// прямой ход метода Гаусса
	//измерение времени прямого хода
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (k = 0; k < rows; ++k)
	{
		//использование cilk_for во внутреннем цикле прямого хода
		for (int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];
			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> duration = (t2 - t1);
	cout << "Consistent realization" << endl;
	cout << "Forward elimination duration for 1500-lines matrix is " << duration.count() << " seconds" << endl;

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];
		//использование cilk_for во внутреннем цикле обратного хода
		for (int j = k + 1; j < rows; ++j)
		{
			result[k] -= matrix[k][j] * result[j];
		}

		result[k] /= matrix[k][k];
	}
	return (t2 - t1);
}

duration<double> SerialGaussMethod_parall(double **matrix, const int rows, double* result)
{
	int k;
	double koef;

	// прямой ход метода Гаусса
	//измерение времени прямого хода
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (k = 0; k < rows; ++k)
	{
		//использование cilk_for во внутреннем цикле прямого хода
		cilk_for(int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];
			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> duration = (t2 - t1);
	cout << "Parallel realization" << endl;
	cout << "Forward elimination duration for 1500-lines matrix is " << duration.count() << " seconds" << endl;

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];
		//использование cilk_for во внутреннем цикле обратного хода
		cilk::reducer_opadd <double> temp(result[k]);
		cilk_for(int j = k + 1; j < rows; ++j)
		{
			temp -= matrix[k][j] * result[j];
		}
		result[k] = temp.get_value();
		result[k] /= matrix[k][k];
	}
	return (t2 - t1);
}


int main()
{
	srand((unsigned)time(0));

	int i;

	// кол-во строк в матрице, приводимой в качестве примера
	//const int test_matrix_lines = 4;

	double **test_matrix = new double*[MATRIX_SIZE];

	// цикл по строкам
	for (i = 0; i < MATRIX_SIZE; ++i)
	{
		// (test_matrix_lines + 1)- количество столбцов в тестовой матрице,
		// последний столбец матрицы отведен под правые части уравнений, входящих в СЛАУ
		test_matrix[i] = new double[MATRIX_SIZE + 1];
	}

	// массив решений СЛАУ
	double *result = new double[MATRIX_SIZE];

	//Заполняем матрицу случайными числами
	InitMatrix(test_matrix);

	//отношение времени выполнения методов
	duration<double> time1 = SerialGaussMethod_posled(test_matrix, MATRIX_SIZE, result);
	duration<double> time2 = SerialGaussMethod_parall(test_matrix, MATRIX_SIZE, result);
	double ratio = time1 / time2;
	cout << "Time ratio is " << ratio << endl;





	for (i = 0; i < MATRIX_SIZE; ++i)
	{
		delete[]test_matrix[i];
	}
	delete[] result;

	return 0;
}