// Copyright 2021 Alibekov Murad

#include <vector>
#include <complex>
#include <iostream>
#include "./ccs_complex_matrix.h"

int main() {
    ccs_matrix tmp;
    InitializeMatrix(4, 6, tmp);
    
    tmp.values[0] = 9;
    tmp.values[1] = 3;
    tmp.values[2] = 8;
    tmp.values[3] = 15;
    tmp.values[4] = 7;
    tmp.values[5] = 16;
    
    tmp.rows[0] = 3;
    tmp.rows[1] = 0;
    tmp.rows[2] = 1;
    tmp.rows[3] = 3;
    tmp.rows[4] = 0;
    tmp.rows[5] = 3;
    
    tmp.col_indexes[0] = 0;
    tmp.col_indexes[1] = 1;
    tmp.col_indexes[2] = 2;
    tmp.col_indexes[3] = 4;
    tmp.col_indexes[4] = 6;
    
    
    PrintDensificationOfCCSMatrix(tmp);
    PrintDensificationOfCCSMatrix(transpose(tmp));
    // PrintMatrixT(4, 6, tmp);
    PrintDensificationOfCCSMatrix(naive_multiplicate(tmp, transpose(tmp)));
    
    
    return 0;
}