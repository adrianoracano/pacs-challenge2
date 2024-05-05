
#include "Matrix.hpp"

namespace algebra {

    template<typename T, StorageOrder Order>
    std::vector<T> operator*(const Matrix<T, Order>& mat, const std::vector<T>& vec) {
    
        std::vector<T> result(mat.rows, 0);

        //if Order == StorageOrder::ROW_MAJOR, the first array index is the row 
        if constexpr (Order == StorageOrder::ROW_MAJOR) {
            if (!mat.compressed) {
                for(auto element : mat.data){
                    auto& i = (element.first)[0]; 
                    auto& j = (element.first)[1];
                    auto& value = element.second;
                    result[i] += value * vec[j]; 
                }
            }else{

                // Perform matrix-vector multiplication for row-major storage
                for (std::size_t i = 0; i < mat.rows; ++i) {
                    for (std::size_t j = mat.row_indices[i]; j < mat.row_indices[i + 1]; ++j) {
                        
                        //j here is the index of the vector of values
                        //the actual j (column) is mat.col_indices[j]], since we are in CSR

                        result[i] += mat.values[j] * vec[mat.col_indices[j]];
                    }
                }
            }
        }else{
            if (!mat.compressed) {
                for(auto element : mat.data){
                    auto& j = (element.first)[0]; //i and j are inverted wrt to the ROW_MAJOR case
                    auto& i = (element.first)[1];
                    auto& value = element.second;
                    result[i] += value * vec[j]; 
                }
            }else{

                //Weighted sum of columns
                for (std::size_t j = 0; j < mat.cols; ++j) {
                    for (std::size_t i = mat.col_indices[j]; i < mat.col_indices[j + 1]; ++i) {

                        //i here is the index of the vector of values
                        //the actual i (row) is mat.row_indices[i], since we are in CSC
                        //this is equivalent to looping along the rows and doing a weighted sum of the columns

                        result[mat.row_indices[i]] += mat.values[i] * vec[j];
                    }
                }
            }
        }
        return result;    
    }   

    // Explicit instantiation for matrix-vector multiplication
    template std::vector<double> operator*<double, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::ROW_MAJOR>&, const std::vector<double>&);
    template std::vector<double> operator*<double, StorageOrder::COLUMN_MAJOR>(const Matrix<double, StorageOrder::COLUMN_MAJOR>&, const std::vector<double>&);
    template std::vector<std::complex<double>> operator*<std::complex<double>, StorageOrder::ROW_MAJOR>(const Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>&, const std::vector<std::complex<double>>&);
    template std::vector<std::complex<double>> operator*<std::complex<double>, StorageOrder::COLUMN_MAJOR>(const Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>&, const std::vector<std::complex<double>>&);

}