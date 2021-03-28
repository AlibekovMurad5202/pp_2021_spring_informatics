// Copyright 2021 Alibekov Murad

/* Это не работающая версия кода, т.е. 
 * по факту этот файл является лишь
 * неким сборником обрывков кода для случая,
 * если придётся по той или иной причине
 * переделывать в C-style код
 */


#include <vector>
#include <complex>
#include <iostream>
#include <cstring>


void PrintMatrixT(int N, int NZ, ccs_matrix &mtx)
{
    for (int i = 0; i < N; i++) {
        for (int j = mtx.col_indexes[i]; j < mtx.col_indexes[i+1]; j++) {
            std::cout << mtx.values[j].real() << "[" << mtx.rows[j] << "," << i << "] ";
            /*for (int k = 0; k < N; k++) {
                if (mtx.rows[k] == k) {
                    std::cout << mtx.values[j].real() << "[" << k << "," << j << "] ";
                } else {
                    std::cout << 0 << "[" << k << "," << j << "] ";
                }
                // std::cout << "[" << mtx.rows[j] << "," << j << "] ";
            }*/
        }
        std::cout << std::endl;
    }
};

void PrintDensificationOfCCSMatrix(ccs_matrix &mtx)
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
            std::cout << dense_matrix[i][j].real() << " ";
        }
        std::cout << std::endl;
    }
};

ccs_matrix naive_multiplicate(const ccs_matrix &A, const ccs_matrix &B)
{
    ccs_matrix AT = transpose(A);
    if (A.N != B.N) 
    {
        // ouch!
    };
    
    const double ZERO_IN_CCS = 0.00001;
    int N = A.N;
    
    std::vector<int> rows;
    std::vector<std::complex<double> > values;
    std::vector<int> col_indexes;
    //clock_t start = clock();

    col_indexes.push_back(0);
    for (int i = 0; i < N; i++)
    {
        int colNZ = 0;
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
                values.push_back(sum);
                colNZ++;
                // std::cout << " " << sum << "[" << j << ", ]";
            }
        }
        col_indexes.push_back(colNZ + col_indexes[i]);
    }
    
    ccs_matrix C;
    InitializeMatrix(N, rows.size(), C);
    for (int j = 0; j < rows.size(); j++) {
        C.rows[j] = rows[j];
        C.values[j] = values[j];
    }
    for(int i = 0; i <= N; i++)
        C.col_indexes[i] = col_indexes[i];
    //clock_t finish = clock();
    //time = (double)(finish - start) / CLOCKS_PER_SEC;
    return C;
};



struct ccs_matrix {
    int N;
    int NZ;

    std::complex<double>* values;
    int* rows;
    int* col_indexes;
}

void InitializeMatrix(int N, int NZ, ccs_matrix &mtx)
{
    mtx.N = N;
    mtx.NZ = NZ;
    mtx.values = new std::complex<double>[NZ];
    mtx.rows = new int[NZ];
    mtx.col_indexes = new int[N + 1];
} 

ccs_matrix transpose(const css_matrix &A) {
    ccs_matrix AT;
    InitializeMatrix(A.N, A.NZ, AT);
    
    memset(AT.col_indexes, 0, (N+1) * sizeof(int));
    for (i = 0; i < A.NZ; i++)
        AT.col_indexes[A.rows[i] + 1]++;
    
    int S = 0;
    for (i = 1; i <= A.N; i++)
    {
        tmp = AT.col_indexes[i];
        AT.col_indexes[i] = S;
        S = S + tmp;
    }
    
    for (i = 0; i < A.N; i++)
    {
        int j1 = A.col_indexes[i]; j2 = A.col_indexes[i+1];
        AT_row = i; // Строка в AT - столбец в А
        for (int j = j1; j < j2; j++)
        {
            std::complex<double> AT_V = A.values[j]; // Значение
            int AT_col_index = A.rows[j]; // Столбец в AT
            int AT_j_index = AT.col_indexes[AT_col_index + 1];
            AT.values[AT_j_index] = AT_V;
            AT.rows[AT_j_index] = AT_row;
            AT.col_indexes[AT_col_index + 1]++;
        }
    }
    
    return AT;
}

void FreeMatrix(ccs_matrix &mtx)
{
 delete[] mtx.values;
 delete[] mtx.rows;
 delete[] mtx.col_indexes;
}

#endif  // MODULES_TASK_1_ALIBEKOV_MURAD_CCS_COMPLEX_MATRIX_CCS_COMPLEX_MATRIX_H_