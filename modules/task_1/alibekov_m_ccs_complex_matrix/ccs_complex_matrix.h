// Copyright 2021 Alibekov Murad
#ifndef MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_
#define MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_

#include <vector>
#include <complex>
#include <iostream>

#include <cstdlib> //srand


const double ZERO_IN_CCS = 0.00001;
const double MAX_VAL = 1000;

struct ccs_matrix {
    int N;       // size of matrix (N x N)
    int NZ;      // count of non-zero elements
    
    // array of values (size = NZ):
    std::vector<std::complex<double> > values;
    
    // array of rows' numbers (size = NZ):
    std::vector<int> rows;
    
    // array of columns' indexes (size = N + 1):
    std::vector<int> col_indexes;
    
    ccs_matrix(int _N, int _NZ) {
        N = _N;
        NZ = _NZ;
        values = std::vector<std::complex<double> >(NZ);
        rows = std::vector<int>(NZ);
        col_indexes = std::vector<int>(N + 1);
    }
};

void PrintCCSMatrix(ccs_matrix &mtx, bool isComplex=false)
{
    std::cout << "Matrix [" << &mtx << "] : \n\tvalues: [ ";
    for (int i = 0; i < mtx.NZ; i++)
        std::cout << (isComplex ? mtx.values[i] : mtx.values[i].real()) << " ";
    std::cout << "]\n\trows: [ ";
    for (int i = 0; i < mtx.NZ; i++)
        std::cout << mtx.rows[i] << " ";
    std::cout << "]\n\tcol_indexes: [ ";
    for (int i = 0; i < mtx.N + 1; i++)
        std::cout << mtx.col_indexes[i] << " ";
    std::cout << "]\n";
};

void PrintDensificationOfCCSMatrix(ccs_matrix &mtx, bool isComplex=false)
{
    int N = mtx.N;
    std::vector<std::vector<std::complex<double> > > dense_matrix(N);
    for (int i = 0; i < N; i++)
        dense_matrix[i] = std::vector<std::complex<double> >(N);
    for (int i = 0; i < N; i++)
        for (int j = mtx.col_indexes[i]; j < mtx.col_indexes[i+1]; j++)
            dense_matrix[mtx.rows[j]][i] = mtx.values[j];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            isComplex ? std::cout << dense_matrix[i][j] << " " :
            std::cout << dense_matrix[i][j].real() << " ";
        }
        std::cout << std::endl;
    }
};

ccs_matrix transpose(const ccs_matrix &A) {
    ccs_matrix AT(A.N, A.NZ);
    
    for (int i = 0; i < A.NZ; i++)
        AT.col_indexes[A.rows[i] + 1]++;
    
    int S = 0;
    for (int i = 1; i <= A.N; i++)
    {
        int tmp = AT.col_indexes[i];
        AT.col_indexes[i] = S;
        S = S + tmp;
    }
    
    for (int i = 0; i < A.N; i++)
    {
        int AT_row = i; // Строка в AT - столбец в А
        for (int j = A.col_indexes[i]; j < A.col_indexes[i+1]; j++)
        {
            std::complex<double> AT_V = A.values[j]; // Значение
            int AT_col_index = A.rows[j]; // Столбец в AT - строка в A
            int AT_j_index = AT.col_indexes[AT_col_index + 1];
            AT.values[AT_j_index] = AT_V;
            AT.rows[AT_j_index] = AT_row;
            AT.col_indexes[AT_col_index + 1]++;
        }
    }
    
    return AT;
};

ccs_matrix naive_multiplicate(const ccs_matrix &A, const ccs_matrix &B)
{
    ccs_matrix AT = transpose(A);
    if (A.N != B.N) 
    {
        // ouch!
    };
    
    int N = A.N;
    
    int rows_count = 0;
    std::vector<int> rows;
    std::vector<std::complex<double> > values;
    std::vector<int> col_indexes;

    col_indexes.push_back(0);
    for (int i = 0; i < N; i++)
    {
        int count_NZ = 0;
        for (int j = 0; j < N; j++)
        {
            std::complex<double> sum = {0, 0};
            
            // dot_product
            for (int k = AT.col_indexes[i]; k < AT.col_indexes[i + 1]; k++)
                for (int l = B.col_indexes[j]; l < B.col_indexes[j + 1]; l++)
                    if (AT.rows[k] == B.rows[l])
                    {
                        sum += AT.values[k] * B.values[l];
                        break;
                    } 
            
            
            if ((fabs(sum.real()) > ZERO_IN_CCS) || (fabs(sum.imag()) > ZERO_IN_CCS)) {
                rows.push_back(j);
                rows_count++;
                values.push_back(sum);
                count_NZ++;
            }
        }
        col_indexes.push_back(count_NZ + col_indexes[i]);
    }
    
    ccs_matrix C(N, rows_count);
    for (int j = 0; j < rows.size(); j++) {
        C.rows[j] = rows[j];
        C.values[j] = values[j];
    }
    for(int i = 0; i <= N; i++)
        C.col_indexes[i] = col_indexes[i];
    return C;
};

// Флаг - был ли инициализирован генератор случайных чисел
bool isSrandCalled = false;

double next() { return ((double)rand() / (double)RAND_MAX); };

void swap_int(int &a, int &b) { int tmp = a; a = b; b = tmp; };

// Генерирует квадратную матрицу в формате CCS
// (3 массива, индексация с нуля)
// В каждом столбце count_in_col ненулевых элементов
ccs_matrix generate_regular_ccs(int seed, int N, int count_in_col)
{
    ccs_matrix random_matrix(N, count_in_col * N);
    if (!isSrandCalled)
    {
        srand(seed);
        isSrandCalled = true;
    }
    
    for(int i = 0; i < N; i++)
    {
        // Формируем номера столбцов в столбце i
        for(int j = i * count_in_col; j < (i + 1) * count_in_col; j++)
        {
            // флаг isFound нужен, чтобы мы не наткунулись на ситуацию, 
            // когда 2 элемента рандомно попали на одно и тоже место в матрице
            bool isFound = false;
            do
            {
                random_matrix.rows[j] = rand() % N;
                isFound = true;
                for (int k = i * count_in_col; k < j; k++)
                    if (random_matrix.rows[j] == random_matrix.rows[k])
                        isFound = false;
            } while (!isFound);
        }
        // Сортируем (BubbleSort) номера строк в столбце i
        for (int j = 0; j < count_in_col - 1; j++)
            for (int k = i * count_in_col; k < (i + 1) * count_in_col - 1 - j; k++)
                if (random_matrix.rows[k] > random_matrix.rows[k + 1])
                    swap_int(random_matrix.rows[k], random_matrix.rows[k + 1]);
    }
    
    // Заполняем массив значений
    for (int i = 0; i < count_in_col * N; i++)
        random_matrix.values[i] = { next() * MAX_VAL, next() * MAX_VAL };
        
    // Заполняем массив индексов столбцов
    for (int i = 0; i <= N; i++)
        random_matrix.col_indexes[i] = i * count_in_col;
    
    return random_matrix;
}










/*
// Генерирует квадратную матрицу в формате CCS
// (3 массива, индексация с нуля)
// Число ненулевых элементов в строках растет
// от 1 до max_count_in_col. Закон роста - кубическая парабола
ccs_matrix generate_special_ccs(int seed, int N, int max_count_in_col)
{
    if (!isSrandCalled)
    {
        srand(seed);
        isSrandCalled = true;
    }
    
    double end_count = pow((double)max_count_in_col, 1.0 / 3.0);
    double step = end_count / N;
    std::vector<std::vector<int> > columns(N);
    int NZ = 0;
    for (int i = 0; i < N; i++)
    {
        int count_in_col = int(pow((double(i + 1) * step), 3) + 1);
        NZ += count_in_col;
        int num1 = (count_in_col - 1) / 2;
        int num2 = count_in_col - 1 - num1;
        if (count_in_col != 0)
        {
            if (i < num1)
            {
                num2 += num1 - i;
                num1 = i;
                for(int j = 0; j < i; j++)
                    columns[i].push_back(j);
                columns[i].push_back(i);
                for(int j = 0; j < num2; j++)
                    columns[i].push_back(i + 1 + j);
            }
            else
            {
                if (N - i - 1 < num2)
                {
                    num1 += num2 - (N - 1 - i);
                    num2 = N - i - 1;
                }
                for (int j = 0; j < num1; j++)
                    columns[i].push_back(i - num1 + j);
                columns[i].push_back(i);
                for (int j = 0; j < num2; j++)
                    columns[i].push_back(i + j + 1);
            }
        }
    }
    
    ccs_matrix random_matrix(N, NZ);
    int count = 0;
    int sum = 0;
    for (int i = 0; i < N; i++)
    {
        random_matrix.col_indexes[i] = sum;
        sum += columns[i].size();
        for (unsigned int j = 0; j < columns[i].size(); j++)
        {
            random_matrix.rows[count] = columns[i][j];
            random_matrix.values[count] = next();
            count++;
        }
    }
    random_matrix.col_indexes[N] = sum;
    
    return random_matrix;
}*/

#endif  // MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_