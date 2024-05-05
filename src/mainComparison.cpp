#include "Matrix.hpp"
#include "MatrixMarketReader.hpp"
#include <chrono>
#include <iostream>
#include <iomanip> // Include <iomanip> for formatting
#include <random>

using namespace algebra;

int main() {

    //Comparison of result and time taken for different state and storage order. 

    // Read matrix from file
    std::string filename = "lnsp_131.mtx";
    Matrix<double, StorageOrder::ROW_MAJOR> matrix_row = readMatrixMarket<double, StorageOrder::ROW_MAJOR>(filename);
    Matrix<double, StorageOrder::COLUMN_MAJOR> matrix_col = readMatrixMarket<double, StorageOrder::COLUMN_MAJOR>(filename);

    // Generate vector of the right dimension
    std::vector<double> vec(matrix_row.getCols(), 1.0); // Fill with ones

    // Time the product in uncompressed state
    auto start_uncompressed_row = std::chrono::high_resolution_clock::now();
    auto result_uncompressed_row = matrix_row * vec;
    auto end_uncompressed_row = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_uncompressed_row = end_uncompressed_row - start_uncompressed_row;

    auto start_uncompressed_col = std::chrono::high_resolution_clock::now();
    auto result_uncompressed_col = matrix_col * vec;
    auto end_uncompressed_col = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_uncompressed_col = end_uncompressed_col - start_uncompressed_col;

    // Compress the matrix
    matrix_row.compress();
    matrix_col.compress();

    // Time the product in compressed state
    auto start_compressed_row = std::chrono::high_resolution_clock::now();
    auto result_compressed_row = matrix_row * vec; // Assuming the matrix is already compressed
    auto end_compressed_row = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_compressed_row = end_compressed_row - start_compressed_row;

    auto start_compressed_col = std::chrono::high_resolution_clock::now();
    auto result_compressed_col = matrix_col * vec; // Assuming the matrix is already compressed
    auto end_compressed_col = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_compressed_col = end_compressed_col - start_compressed_col;

    // Compare the results of uncompressed and compressed products
    std::cout << std::setw(30) << std::right << "ROW MAJOR" << std::setw(40) << std::right << "COLUMN MAJOR" << std::endl;
    std::cout << std::setw(10) << std::left << "Index" << std::setw(20) << std::left << "Uncompressed" 
                                                       << std::setw(20) << std::left << "Compressed"
                                                       << std::setw(20) << std::left << "Uncompressed" 
                                                       << std::setw(20) << std::left << "Compressed"
                                                       << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    for (size_t i = 0; i < result_uncompressed_row.size(); ++i) {
        std::cout << std::setw(20) << std::left << i << std::setw(20) << std::left << result_uncompressed_row[i] 
                                                     << std::setw(20) << std::left << result_compressed_row[i] 
                                                     << std::setw(20) << std::left << result_uncompressed_col[i] 
                                                     << std::setw(20) << std::left << result_compressed_col[i]
                                                     << std::endl;
    }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::setw(20) << std::left << "Time taken" << std::setw(20) << std::left << duration_uncompressed_row.count() 
                                                       << std::setw(20) << std::left << duration_compressed_row.count()
                                                       << std::setw(20) << std::left << duration_uncompressed_col.count() 
                                                       << std::setw(20) << std::left << duration_compressed_col.count()
                                                       << std::endl;

    return 0;
}

