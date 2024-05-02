#include "Matrix.hpp"

namespace algebra {

    template<typename T, StorageOrder Order>
    std::vector<T> operator*(const Matrix<T, Order>& mat, const std::vector<T>& vec) {
    
        std::vector<T> result(mat.rows, 0);

        //if Order == StorageOrder::ROW_MAJOR, the first array index is the row while the second is the column
        if constexpr (Order == StorageOrder::ROW_MAJOR) {
            if (!mat.compressed) {
                for(auto element : mat.data){
                    auto& i = (element.first)[0]; 
                    auto& j = (element.first)[1];
                    auto& value = element.second;
                    result[i] += value * vec[j]; 
                    //std::cout<<result[i]<<" ";
                }
            }else{
                // Perform matrix-vector multiplication for row-major storage (CSR) (my version)
                for (std::size_t i = 0; i < mat.rows; ++i) {
                    for (std::size_t j = mat.row_indices[i]; j < mat.row_indices[i + 1]; ++j) {
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
            // Perform matrix-vector multiplication for column-major storage
                /* for (std::size_t j = 0; j < mat.cols; ++j) {
                    for (std::size_t i = mat.col_indices[j]; i < mat.col_indices[j + 1]; ++i) {
                        result[j] += mat.values[i] * vec[mat.row_indices[i]];
                    }
                } */

                //Weighted sum of columns
                for (std::size_t j = 0; j < mat.cols; ++j) {
                    for (std::size_t i = mat.col_indices[j]; i < mat.col_indices[j + 1]; ++i) {
                        //result[i] += mat.values[] * vec[mat.row_indices[i]];
                        result[mat.row_indices[i]] += mat.values[i] * vec[j];
                    }
                }
            }
        }
        return result;    
    }   

    //Product with 1-column matrices
    /* template<typename T, StorageOrder Order1, StorageOrder Order2>
    std::vector<T> operator*(const Matrix<T, Order1>& mat1, const Matrix<T, Order2>& mat2) {

        if (mat2.cols > 1) 
            throw std::out_of_range("Only 1-column matrices are accepted");

        if (mat2.rows != mat1.cols) 
            throw std::out_of_range("Incompatible sizes");    

        std::vector<T> result(mat1.rows, 0);

        if constexpr (Order1 == StorageOrder::ROW_MAJOR && Order2 == StorageOrder::ROW_MAJOR){
            if(mat1.compressed && mat2.compressed){
                for(auto el1 : mat1.data){
                    auto i1 = (el1.first)[0];
                    auto j1 = (el1.first)[1];
                    for (auto el2 : mat2.data){
                        auto i2 = (el2.first)[0];
                        auto j2 = (el2.first)[1];
                        if (j1 == i2){
                            result[i1] += el1.second * el2.second;
                        }
                    }
                }
            }

            if(!mat1.compressed && !mat2.compressed){


                
                std::size_t i1 = 0;
                for(std::size_t val_id1 = 0 ; val_id1 < mat1.values.size(); val_id1++){
                    double val1 = mat1.values[val_id1];
                    while(mat1.row_indices[i1+1] < val_id1 + 1){
                        i1++; 
                    }
                    std::size_t j1 = mat1.col_indices[val1];
                    for (auto el2 : mat2.data){
                        auto i2 = (el2.first)[0];
                        auto j2 = (el2.first)[1];
                        if (j1 == i2){
                            result[i1] += el1.second * el2.second;
                        }
                    }
                    count++;
                }
            }
        }

        return result;    
    } */

    // Explicit instantiation for matrix-vector multiplication
    template std::vector<double> operator*<double, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::ROW_MAJOR>&, const std::vector<double>&);
    template std::vector<double> operator*<double, StorageOrder::COLUMN_MAJOR>(const Matrix<double, StorageOrder::COLUMN_MAJOR>&, const std::vector<double>&);

}