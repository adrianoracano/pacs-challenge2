#include "Matrix.hpp"
#include <fstream>
#include <sstream>

namespace algebra{
    template<typename T, StorageOrder Order>
    Matrix<T, Order> readMatrixMarket(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        std::string line;
        std::size_t rows = 0, cols = 0, nnz = 0;
        bool header_read = false;

        // Read the header to get matrix dimensions and number of non-zero entries
        while (std::getline(file, line) && !header_read) {
            if (line[0] != '%') {
                std::istringstream iss(line);
                iss >> rows >> cols >> nnz;
                header_read = true;
            }
        }

        if (!header_read) {
            throw std::runtime_error("Failed to read Matrix Market header from file: " + filename);
        }

        // Create a matrix with the read dimensions
        Matrix<T, Order> matrix(rows, cols); //by default in uncompressed state

        // Read the non-zero entries and insert them into the matrix
        for (std::size_t i = 0; i < nnz; ++i) {
            std::size_t row, col;
            T value;
            file >> row >> col >> value;
            matrix(row - 1, col - 1) = value; // Matrix Market indices are 1-based
        }

        return matrix;
    }
}