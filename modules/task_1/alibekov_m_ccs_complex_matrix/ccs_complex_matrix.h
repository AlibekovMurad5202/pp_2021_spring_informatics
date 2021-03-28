// Copyright 2021 Alibekov Murad
#ifndef MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_
#define MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_

#include <vector>
#include <complex>
#include <iostream>

const double ZERO_IN_CCS = 0.00001;
const double MAX_VAL = 1000;

static bool isSrandCalled = false;

struct ccs_complex_matrix {
    int N;       // size of matrix (N x N)
    int NZ;      // count of non-zero elements
    
    // array of values (size = NZ):
    std::vector<std::complex<double> > values;
    
    // array of rows' numbers (size = NZ):
    std::vector<int> rows;
    
    // array of columns' indexes (size = N + 1):
    std::vector<int> col_indexes;
    
    ccs_complex_matrix(int _N, int _NZ) {
        N = _N;
        NZ = _NZ;
        values = std::vector<std::complex<double> >(NZ);
        rows = std::vector<int>(NZ);
        col_indexes = std::vector<int>(N + 1);
    }
    
    friend bool operator==(const ccs_complex_matrix &A, const ccs_complex_matrix &B);
};

ccs_complex_matrix generate_regular_ccs(int seed, int N, int count_in_col);

ccs_complex_matrix transpose(const ccs_complex_matrix &A);
ccs_complex_matrix naive_multiplicate(const ccs_complex_matrix &A, const ccs_complex_matrix &B);
ccs_complex_matrix optim_multiplicate(const ccs_complex_matrix &A, const ccs_complex_matrix &B);

void PrintCCSMatrix(const ccs_complex_matrix &A, bool isComplex=true);
void PrintDensificationOfCCSMatrix(const ccs_complex_matrix &A, bool isComplex=true);

inline double next() { return ((double)rand() / (double)RAND_MAX); };
inline void swap_int(int &a, int &b) { int tmp = a; a = b; b = tmp; };

bool operator==(const ccs_complex_matrix &A, const ccs_complex_matrix &B);

#endif  // MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_