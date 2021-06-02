// Copyright 2021 Alibekov Murad
#include <gtest/gtest.h>
#include <vector>
#include <complex>
#include <iostream>
#include <ctime>
#include "./ccs_complex_matrix.h"

int SEED_1 = 86538;
int SEED_2 = 2395;
int N = 5000;
int COUNT_IN_COL = 100;
int TBB_GRANSIZE = 1;


TEST(SPARSE_MATRICES, PRINT_SPARSE_MATRIX) {
    ccs_complex_matrix sparse_matrix(4, 6);
    sparse_matrix.values = { 9, 3, 8, 15, 7, 16 };
    sparse_matrix.rows = { 3, 0, 1, 3, 0, 3 };
    sparse_matrix.col_indexes = { 0, 1, 2, 4, 6 };

    std::cout << "\tFirst type:\n";
    PrintCCSMatrix(sparse_matrix);
    std::cout << std::endl;

    std::cout << "\tSecond type:\n";
    PrintCCSMatrix(sparse_matrix, false);
    std::cout << std::endl;

    std::cout << "\tThird type:\n";
    PrintDensificationOfCCSMatrix(sparse_matrix);
    std::cout << std::endl;

    std::cout << "\tFourth type:\n";
    PrintDensificationOfCCSMatrix(sparse_matrix, false);
    std::cout << std::endl;
}


/////////////////////////////////////////////
///    NAIVE_MULTIPLY_SPARSE_MATRICES    ////
/////////////////////////////////////////////

TEST(NAIVE_MULTIPLY_SPARSE_MATRICES, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    tbb::tick_count start_time, finish_time;
    start_time = tbb::tick_count::now();
    EXPECT_NO_THROW(naive_multiplicate(big_sparse_matrix_1, big_sparse_matrix_2));
    finish_time = tbb::tick_count::now();

    printf("\tTime  = %f\n", (finish_time - start_time).seconds());
}


/////////////////////////////////////////////
///    OPTIM_MULTIPLY_SPARSE_MATRICES    ////
/////////////////////////////////////////////

TEST(OPTIM_MULTIPLY_SPARSE_MATRICES, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    tbb::tick_count start_time, finish_time;
    start_time = tbb::tick_count::now();
    EXPECT_NO_THROW(optim_multiplicate(big_sparse_matrix_1, big_sparse_matrix_2));
    finish_time = tbb::tick_count::now();

    printf("\tTime  = %f\n", (finish_time - start_time).seconds());
}


/////////////////////////////////////////////////
///    NAIVE_MULTIPLY_SPARSE_MATRICES_OMP    ////
/////////////////////////////////////////////////

TEST(NAIVE_MULTIPLY_SPARSE_MATRICES_OMP, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    double start_time = omp_get_wtime();
    EXPECT_NO_THROW(naive_multiplicate_omp(big_sparse_matrix_1, big_sparse_matrix_2));
    double finish_time = omp_get_wtime();

    printf("\tTime  = %f\n", finish_time - start_time);
}


/////////////////////////////////////////////////
///    OPTIM_MULTIPLY_SPARSE_MATRICES_OMP    ////
/////////////////////////////////////////////////

TEST(OPTIM_MULTIPLY_SPARSE_MATRICES_OMP, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    double start_time = omp_get_wtime();
    EXPECT_NO_THROW(optim_multiplicate_omp(big_sparse_matrix_1, big_sparse_matrix_2));
    double finish_time = omp_get_wtime();

    printf("\tTime  = %f\n", finish_time - start_time);
}


/////////////////////////////////////////////////
///    NAIVE_MULTIPLY_SPARSE_MATRICES_TBB    ////
/////////////////////////////////////////////////

TEST(NAIVE_MULTIPLY_SPARSE_MATRICES_TBB, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    tbb::tick_count start_time, finish_time;
    start_time = tbb::tick_count::now();
    EXPECT_NO_THROW(naive_multiplicate_tbb(big_sparse_matrix_1, big_sparse_matrix_2, TBB_GRANSIZE));
    finish_time = tbb::tick_count::now();

    printf("\tTime  = %f\n", (finish_time - start_time).seconds());
}


/////////////////////////////////////////////////
///    OPTIM_MULTIPLY_SPARSE_MATRICES_TBB    ////
/////////////////////////////////////////////////

TEST(OPTIM_MULTIPLY_SPARSE_MATRICES_TBB, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    tbb::tick_count start_time, finish_time;
    start_time = tbb::tick_count::now();
    EXPECT_NO_THROW(optim_multiplicate_tbb(big_sparse_matrix_1, big_sparse_matrix_2, TBB_GRANSIZE));
    finish_time = tbb::tick_count::now();

    printf("\tTime  = %f\n", (finish_time - start_time).seconds());
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
