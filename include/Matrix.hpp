
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <map>
#include <array>
#include <string>
#include <cstddef>
#include <stdexcept>
#include <iostream>

namespace algebra {

    enum class StorageOrder { ROW_MAJOR, COLUMN_MAJOR };
    enum class NormType {One, Frobenius, Infinity};

    template<typename T, StorageOrder Order>
    class Matrix {
    public:
        Matrix();
        explicit Matrix(std::size_t size);
        Matrix(std::size_t rows, std::size_t cols);

        void resize(std::size_t rows, std::size_t cols);
        void compress();
        void uncompress();
        bool isCompressed() const;

        T operator()(std::size_t i, std::size_t j) const;
        T& operator()(std::size_t i, std::size_t j);

        template<typename U, StorageOrder S>
        friend std::vector<U> operator*(const Matrix<U, S>& mat, const std::vector<U>& vec);
        
        template<typename U, StorageOrder Order1, StorageOrder Order2>
        friend std::vector<U> operator*(const Matrix<U, Order1>& mat1, const Matrix<U, Order2>& mat2);

        template<typename U, StorageOrder S, NormType N>
        friend U norm(const Matrix<U, S>& mat);
        //friend T norm<T, Order, NormType::Frobenius>(const Matrix<T, Order>& mat);
        //friend T norm<T, Order, NormType::Infinity>(const Matrix<T, Order>& mat);
        
        template<NormType N>
        double norm();

        //getters
        std::size_t getRows() const {
            return rows;
        }

        std::size_t getCols() const {
            return cols;
        }

        void PrintRowIndices(){
            if(!compressed)
                throw std::runtime_error("Not in compressed format.");
            else{
                std::cout<<"CSR Row indices (size = "<< row_indices.size() <<"):"<<std::endl;
                for(std::size_t i = 0; i < row_indices.size(); ++i){
                    std::cout<<row_indices[i]<<" ";
                }
                std::cout<<std::endl;
            }
        }

        void PrintColumnIndices(){
            if(!compressed)
                throw std::runtime_error("Not in compressed format.");
            else{
                std::cout<<"CSR Column indices (size = "<< col_indices.size() << "):"<<std::endl;
                for(std::size_t i = 0; i < col_indices.size(); ++i){
                    std::cout<<col_indices[i]<<" ";
                }
                std::cout<<std::endl;
            }
        }

        void PrintValues(){
            if(!compressed){
                std::cout<<"COO values (size = " << data.size() << "):"<<std::endl;
                for(auto val : data){
                    std::cout<< val.second <<" ";
                }
                std::cout<<std::endl;
            }else{
                std::cout<<"CSR values (size = " << values.size() << "):"<<std::endl;
                for(std::size_t i = 0; i < values.size(); ++i){
                    std::cout<<values[i]<<" ";
                }
                std::cout<<std::endl;
            }
        }
        

    private:
        std::vector<std::size_t> row_indices; // Row indices for CSR format
        std::vector<std::size_t> col_indices; // Column indices for CSR format
        std::vector<T> values; // Values for CSR format
        std::map<std::array<std::size_t, 2>, T> data; // COOmap format
        bool compressed;
        std::size_t rows;
        std::size_t cols;
    };

    /* template<typename T, StorageOrder Order, NormType Norm>
    T norm(const Matrix<T, Order>& mat); */
    

} // namespace algebra

#endif // MATRIX_HPP
