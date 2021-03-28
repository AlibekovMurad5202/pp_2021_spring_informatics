// Copyright 2021 Alibekov Murad
#ifndef MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_
#define MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_

#include <vector>
#include <complex>
#include <iostream>

const double ZERO_IN_CCS = 0.00001;
const double MAX_VAL = 1000;

// Флаг - был ли инициализирован генератор случайных чисел
static bool isSrandCalled = false;

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

// Генерирует квадратную матрицу в формате CCS
// (3 массива, индексация с нуля)
// В каждом столбце count_in_col ненулевых элементов
ccs_matrix generate_regular_ccs(int seed, int N, int count_in_col);

ccs_matrix transpose(const ccs_matrix &A);
ccs_matrix naive_multiplicate(const ccs_matrix &A, const ccs_matrix &B);

void PrintCCSMatrix(const ccs_matrix &A, bool isComplex=true);
void PrintDensificationOfCCSMatrix(const ccs_matrix &A, bool isComplex=true);

inline double next() { return ((double)rand() / (double)RAND_MAX); };
inline void swap_int(int &a, int &b) { int tmp = a; a = b; b = tmp; };

#endif  // MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_