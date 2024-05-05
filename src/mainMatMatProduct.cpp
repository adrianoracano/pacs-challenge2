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
    std::string filename = "/lnsp_131.mtx";
    Matrix<double, Order> matrix = readMatrixMarket<double, Order>(filename);

    // Print a portion of the matrix
    matrix.print(10,10);
    //printMatrix(matrix, 10, 10); // Print first 5 rows and 5 columns
    std::cout<<"1-Norm of matrix = "<< matrix.norm< NormType::One>()<<std::endl;
    std::cout<<"Infinity-Norm of matrix = "<< matrix.norm< NormType::Infinity>()<<std::endl;
    std::cout<<"Frobenius-Norm of matrix  = "<< matrix.norm< NormType::Frobenius>()<<std::endl;

    //Trying product between 2 matrices
    std::cout << "Product between 2 matrices: " << std::endl;

    std::size_t col2 = 5;
    Matrix<double, StorageOrder::ROW_MAJOR> matrix2(matrix.getRows(), col2);

    for(std::size_t i = 0; i < matrix.getRows() ; ++i){
        for(std::size_t j = 0 ; j < col2 ; j++)
            matrix2(i,j) = 1 * (j+1);
    };

    Matrix<double, StorageOrder::ROW_MAJOR> result = matrix * matrix2;

    for (std::size_t i = 0; i < result.getRows(); ++i) {
        std::cout << std::setw(10) << std::left << i << std::setw(20) << std::left;
        for (std::size_t j = 0; j < result.getCols(); ++j){ 
            std::cout << result(i, j) << std::setw(20) <<std::left;
        }
        std::cout<<std::endl;  
    }
    
    return 0;
}



