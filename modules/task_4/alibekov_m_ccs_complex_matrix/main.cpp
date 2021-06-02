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

    auto start_time = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(naive_multiplicate(big_sparse_matrix_1, big_sparse_matrix_2));
    auto finish_time = std::chrono::high_resolution_clock::now();

    printf("\tTime  = %f\n", static_cast<float>(
        std::chrono::duration_cast<std::chrono::milliseconds>
            (finish_time - start_time).count() / 1000.));
}


/////////////////////////////////////////////
///    OPTIM_MULTIPLY_SPARSE_MATRICES    ////
/////////////////////////////////////////////

TEST(OPTIM_MULTIPLY_SPARSE_MATRICES, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    auto start_time = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(optim_multiplicate(big_sparse_matrix_1, big_sparse_matrix_2));
    auto finish_time = std::chrono::high_resolution_clock::now();

    printf("\tTime  = %f\n", static_cast<float>(
        std::chrono::duration_cast<std::chrono::milliseconds>
            (finish_time - start_time).count() / 1000.));
}


/////////////////////////////////////////////////
///    NAIVE_MULTIPLY_SPARSE_MATRICES_STD    ////
/////////////////////////////////////////////////

TEST(NAIVE_MULTIPLY_SPARSE_MATRICES_STD, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    auto start_time = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(naive_multiplicate_std(big_sparse_matrix_1, big_sparse_matrix_2));
    auto finish_time = std::chrono::high_resolution_clock::now();

    printf("\tTime  = %f\n", static_cast<float>(
        std::chrono::duration_cast<std::chrono::milliseconds>
            (finish_time - start_time).count() / 1000.));
}


/////////////////////////////////////////////////
///    OPTIM_MULTIPLY_SPARSE_MATRICES_STD    ////
/////////////////////////////////////////////////

TEST(OPTIM_MULTIPLY_SPARSE_MATRICES_STD, PERFORMANCE_MEASUREMENT_OF_MULTIPLICATION_BIG_SPARSE_MATRICES) {
    ccs_complex_matrix big_sparse_matrix_1 = generate_regular_ccs(SEED_1, N, COUNT_IN_COL);
    std::cout << "\tFirst matrix is generated!\n";

    ccs_complex_matrix big_sparse_matrix_2 = generate_regular_ccs(SEED_2, N, COUNT_IN_COL);
    std::cout << "\tSecond matrix is generated!\n";

    auto start_time = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(optim_multiplicate_std(big_sparse_matrix_1, big_sparse_matrix_2));
    auto finish_time = std::chrono::high_resolution_clock::now();

    printf("\tTime  = %f\n", static_cast<float>(
        std::chrono::duration_cast<std::chrono::milliseconds>
            (finish_time - start_time).count() / 1000.));
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
