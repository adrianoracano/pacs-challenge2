
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <map>
#include <array>
#include <string>
#include <cstddef>
#include <stdexcept>
#include <iostream>
#include <complex>
#include <iomanip>

namespace algebra {

    
    enum class StorageOrder { ROW_MAJOR, COLUMN_MAJOR }; 
    enum class NormType {One, Frobenius, Infinity};      

    template<typename T, StorageOrder Order>
    class Matrix {
    public:

        //Constructors
        Matrix();
        explicit Matrix(std::size_t size);
        Matrix(std::size_t rows, std::size_t cols);

        //Methods for switching format
        void resize(std::size_t rows, std::size_t cols);
        void compress();
        void uncompress();
        bool isCompressed() const;

        //Const and non-const version of the call operator
        T operator()(std::size_t i, std::size_t j) const;
        T& operator()(std::size_t i, std::size_t j);

        //Matrix-vector product
        template<typename U, StorageOrder S>
        friend std::vector<U> operator*(const Matrix<U, S>& mat, const std::vector<U>& vec);
        
        //Matrix-matrix product
        template<typename U, StorageOrder Order1, StorageOrder Order2>
        friend Matrix<U, Order1> operator*(const Matrix<U, Order1>& mat1, const Matrix<U, Order2>& mat2);

        /* 
        template<typename U, StorageOrder Order1, StorageOrder Order2>
        friend std::vector<U> operator*(const Matrix<U, Order1>& mat1, const Matrix<U, Order2>& mat2); 
        */
        
        //Tempate function for the norm
        template<NormType N>
        double norm();

        //getters
        std::size_t getRows() const {
            return rows;
        }

        std::size_t getCols() const {
            return cols;
        }

        std::size_t getNnz() const {
            if (!compressed)
                return data.size();
            else
                return values.size();
        }

        //Print arbitrary number of rows and cols
        void print(std::size_t n_rows, std::size_t n_cols){
            std::cout << "Matrix dimensions: " << rows << " x " << cols << std::endl;
            std::cout << "Print of the first " << n_rows << " rows and " << n_cols << " columns:" << std::endl;
            std::size_t maxRows = std::min(n_rows, rows);
            std::size_t maxCols = std::min(n_cols, cols);
            for (std::size_t i = 0; i < maxRows; ++i) {
                for (std::size_t j = 0; j < maxCols; ++j) {
                    std::cout << std::setw(10) << std::left << (*this)(i, j) << " ";
                }
                std::cout << std::endl;
            }
        }
        

    private:

        std::vector<std::size_t> row_indices; // Row indices for CSR format
        std::vector<std::size_t> col_indices; // Column indices for CSR format
        std::vector<T> values; // Values for CSR format
        std::map<std::array<std::size_t, 2>, T> data; // COOmap format (in row major the row comes first, in column major it's the opposite)
        bool compressed;
        std::size_t rows;
        std::size_t cols;
    };

} // namespace algebra

#endif // MATRIX_HPP
