#include "Matrix.hpp"
#include "MatrixMarketReader.hpp"
#include <chrono>
#include <iostream>
#include <iomanip> // Include <iomanip> for formatting
#include <random>

using namespace algebra;


int main() {
    constexpr StorageOrder Order = StorageOrder::ROW_MAJOR;

    //Filling randomly the matrix
    std::random_device rd;
    std::mt19937 gen(42); // Mersenne Twister pseudo-random generator
    std::uniform_int_distribution<int> dist(1, 5); // Random integers between 1 and 5

    //Initilizing the matrix and the vector of ones (complex form)
    std::size_t size = 10;
    Matrix<double, Order> matrix_real(size);
    Matrix<std::complex<double>, Order> matrix_complex(size);

    constexpr std::complex<double> one = {1.0, 1.0};
    std::vector<std::complex<double>> v_complex(10, one);

    //Filling with random complex values
    for(std::size_t i = 0; i < 10 ; ++i){
        for(std::size_t j = 0; j < 10 ; ++j){
            matrix_complex(i,j) = {dist(gen), dist(gen)};
        }
    }

    printMatrix(matrix_complex, 10, 10); // Print first 10 rows and 10 columns
    std::cout<<"1-Norm of matrix = "<< matrix_complex.norm< NormType::One>()<<std::endl;
    std::cout<<"Infinity-Norm of matrix = "<< matrix_complex.norm< NormType::Infinity>()<<std::endl;
    std::cout<<"Frobenius-Norm of matrix  = "<< matrix_complex.norm< NormType::Frobenius>()<<std::endl;

    std::vector<std::complex<double>> result_complex = matrix_complex * v_complex;
    
    std::cout<<std::endl;
    std::cout<< " Matrix * v: " << std::endl;
    for (size_t i = 0; i < 10; ++i) {
        std::cout << std::setw(10) << std::left << i << std::setw(20) << std::left << result_complex[i] 
                                                     << std::endl;
    }

    return 0;
}

