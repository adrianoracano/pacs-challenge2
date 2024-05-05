#include "Matrix.hpp"
#include "MatrixMarketReader.hpp"
#include <chrono>
#include <iostream>
#include <iomanip> // Include <iomanip> for formatting
#include <random>

using namespace algebra;

template<typename T, StorageOrder Order>
void printMatrix(const Matrix<T, Order>& matrix, std::size_t rows, std::size_t cols) {
    std::cout << "Matrix dimensions: " << matrix.getRows() << " x " << matrix.getCols() << std::endl;
    std::cout << "Print of the first " << rows << " rows and " << cols << " columns:" << std::endl;
    std::size_t maxRows = std::min(rows, matrix.getRows());
    std::size_t maxCols = std::min(cols, matrix.getCols());
    for (std::size_t i = 0; i < maxRows; ++i) {
        for (std::size_t j = 0; j < maxCols; ++j) {
            std::cout << std::setw(10) << std::left << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }
    
}

int main() {
    // Create an example matrix in uncompressed state
    constexpr StorageOrder Order = StorageOrder::COLUMN_MAJOR;

    // Read matrix from file
    std::string filename = "../test/lnsp_131.mtx";
    Matrix<double, Order> matrix = readMatrixMarket<double, Order>(filename);

    // Print a portion of the matrix
    printMatrix(matrix, 10, 10); // Print first 5 rows and 5 columns
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

    // Output the time taken for the uncompressed product
    std::cout << "Time taken for matrix-vector product (uncompressed): " << duration_uncompressed.count() << " microseconds\n";

    // Compress the matrix (if necessary)
    matrix.compress();

    // Time the product in compressed state
    auto start_compressed = std::chrono::high_resolution_clock::now();
    auto result_compressed = matrix * vec; // Assuming the matrix is already compressed
    auto end_compressed = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_compressed = end_compressed - start_compressed;

    // Output the time taken for the compressed product
    std::cout << "Time taken for matrix-vector product (compressed): " << duration_compressed.count() << " microseconds\n";

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

    //Product with a matrix with just 1 column
    Matrix<double, StorageOrder::COLUMN_MAJOR> v(matrix.getCols(), 1);

    //Filling randomly the matrix
    std::random_device rd;
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_int_distribution<int> dist(1, 5); // Random integers between 1 and 5

    for(std::size_t i = 0; i < v.getRows() ; ++i){
        v(i,0) = 1;
    };

    matrix.uncompress();
    v.compress();

    std::vector<double> result_1col = matrix * v;

    for (size_t i = 0; i < result_1col.size(); ++i) {
        std::cout << std::setw(10) << std::left << i << std::setw(20) << std::left << result_1col[i] 
                                                     << std::endl;

    }

    //Trying product between 2 matrices
    /* std::cout << "Product between 2 matrices: " << std::endl;

    Matrix<double, StorageOrder::ROW_MAJOR> matrix2(matrix.getRows(), 2);

    for(std::size_t i = 0; i < v.getRows() ; ++i){
        for(std::size_t j = 0 ; j < 2 ; j++)
            matrix2(i,j) = 1 * (j+1);
    };

    Matrix<double, StorageOrder::ROW_MAJOR> result_2col = matrix * matrix2;

    for (size_t i = 0; i < result_2col.getRows(); ++i) {
        std::cout << std::setw(10) << std::left << i << std::setw(20) << std::left << result_2col(i, 0) 
                                                                      << std::left << result_2col(i, 1) 
                                                                      << std::endl;

    } */


    

    

    /* //complex case
    constexpr std::complex<double> val = {1.0, 1.0};
    std::vector<std::complex<double>> v_complex(10, val);
    std::size_t size = 10;
    Matrix<double, Order> matrix_real(size);
    Matrix<std::complex<double>, Order> matrix_complex(size);

    //filling with random complex values
    for(std::size_t i = 0; i < 10 ; ++i){
        for(std::size_t j = 0; j < 10 ; ++j){
            matrix_complex(i,j) = {dist(gen), dist(gen)};
        }
    }

    std::cout<<"Complex case: "<<std::endl;

    std::vector<std::complex<double>> result_complex = matrix_complex * v_complex;
    
    for (size_t i = 0; i < 10; ++i) {
        std::cout << std::setw(10) << std::left << i << std::setw(20) << std::left << result_complex[i] 
                                                     << std::endl;
    }

    std::cout<<"1-Norm of complex matrix = "<< matrix_complex.norm< NormType::One>()<<std::endl;
    std::cout<<"Infinity-Norm of complex matrix = "<< matrix_complex.norm< NormType::Infinity>()<<std::endl;
    std::cout<<"Frobenius-Norm of complex matrix  = "<< matrix_complex.norm< NormType::Frobenius>()<<std::endl;
 */
    return 0;
}

/* int main() {

    // Read matrix from file
    std::string filename = "../test/lnsp_131.mtx";
    Matrix<double, StorageOrder::ROW_MAJOR> matrix = readMatrixMarket<double, StorageOrder::ROW_MAJOR>(filename);
    Matrix<double, StorageOrder::COLUMN_MAJOR> matrix2 = readMatrixMarket<double, StorageOrder::COLUMN_MAJOR>(filename);

    matrix2.compress();

    // Generate vector of the right dimension
    std::vector<double> vec(matrix.getCols(), 1.0); // Fill with ones
    auto result = matrix * vec; // Assuming the matrix is already compressed
    auto result2 = matrix2 * vec; // Assuming the matrix is already compressed

    // Compare the results of row-wise and column-wise storage products (if necessary)
    std::cout << std::setw(10) << std::left << "Index" << std::setw(20) << std::left << "Row" 
                                                       << std::setw(20) << std::left << "Column"
                                                       << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        std::cout << std::setw(10) << std::left << i << std::setw(20) << std::left << result[i] 
                                                     << std::setw(20) << std::left << result2[i] 
                                                     << std::endl;
    }


    std::cout<<"1-Norm of matrix 1 = "<< matrix.norm< NormType::One>()<<std::endl;
    std::cout<<"Infinity-Norm of matrix 1 = "<< matrix.norm< NormType::Infinity>()<<std::endl;
    std::cout<<"Frobenius-Norm of matrix 1 = "<< matrix.norm< NormType::Frobenius>()<<std::endl;
    std::cout<<"1-Norm of matrix 2 = "<< matrix2.norm< NormType::One>()<<std::endl;
    std::cout<<"Infinity-Norm of matrix 2 = "<< matrix2.norm< NormType::Infinity>()<<std::endl;
    std::cout<<"Frobenius-Norm of matrix 2 = "<< matrix2.norm< NormType::Frobenius>()<<std::endl;

    return 0;
} */


