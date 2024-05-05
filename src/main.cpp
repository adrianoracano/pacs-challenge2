#include "Matrix.hpp"
#include "MatrixMarketReader.hpp"
#include <chrono>
#include <iostream>
#include <iomanip> // Include <iomanip> for formatting
#include <random>

using namespace algebra;

int main() {
    // Create an example matrix in uncompressed state
    constexpr StorageOrder Order = StorageOrder::ROW_MAJOR;

    // Read matrix from file
    std::string filename = "lnsp_131.mtx";
    Matrix<double, Order> matrix = readMatrixMarket<double, Order>(filename);

    // Print a portion of the matrix
    matrix.print(10,10); // Print first 5 rows and 5 columns
    std::cout<<"1-Norm of matrix = "<< matrix.norm< NormType::One>()<<std::endl;
    std::cout<<"Infinity-Norm of matrix = "<< matrix.norm< NormType::Infinity>()<<std::endl;
    std::cout<<"Frobenius-Norm of matrix  = "<< matrix.norm< NormType::Frobenius>()<<std::endl;

    // Generate vector of the right dimension
    std::vector<double> vec(matrix.getCols(), 1.0); // Fill with ones

    // Time the product in uncompressed state
    auto start_uncompressed = std::chrono::high_resolution_clock::now();
    auto result_uncompressed = matrix * vec;
    auto end_uncompressed = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_uncompressed = end_uncompressed - start_uncompressed;

    // Compress the matrix
    matrix.compress();

    // Time the product in compressed state
    auto start_compressed = std::chrono::high_resolution_clock::now();
    auto result_compressed = matrix * vec; // Assuming the matrix is already compressed
    auto end_compressed = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_compressed = end_compressed - start_compressed;

    // Compare the results of uncompressed and compressed products
    std::cout << std::setw(10) << std::left << "Index" << std::setw(20) << std::left << "Uncompressed" 
                                                       << std::setw(20) << std::left << "Compressed"
                                                       << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    for (size_t i = 0; i < result_uncompressed.size(); ++i) {
        std::cout << std::setw(10) << std::left << i << std::setw(20) << std::left << result_uncompressed[i] 
                                                     << std::setw(20) << std::left << result_compressed[i] 
                                                     << std::endl;
    }

    // Output the time taken for the uncompressed product
    std::cout << "Time taken for matrix-vector product (uncompressed): " << duration_uncompressed.count() << " microseconds\n";
    // Output the time taken for the compressed product
    std::cout << "Time taken for matrix-vector product (compressed): " << duration_compressed.count() << " microseconds\n";
   
    return 0;
}


    