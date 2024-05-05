
#include "Matrix.hpp"

namespace algebra{

    // The case of product with 1-column matrices has been omitted for brevity. It's equivalent 
    // to the below code but the loop along matrix 2 columns is skipped, since only column 0 is 
    // present. Therefore, there's no need to initialize j2 in all cases.

    //Product with general matrices
    template<typename T, StorageOrder Order1, StorageOrder Order2>
    Matrix<T, Order1> operator*(const Matrix<T, Order1>& mat1, const Matrix<T, Order2>& mat2) {

        if (mat2.rows != mat1.cols) 
            throw std::out_of_range("Incompatible sizes");    

        Matrix<T, Order1> result(mat1.rows, mat2.cols);

        if constexpr (Order1 == StorageOrder::ROW_MAJOR && Order2 == StorageOrder::ROW_MAJOR){

            if(!mat1.compressed && !mat2.compressed){
                for(auto el1 : mat1.data){
                    auto& i1 = (el1.first)[0];
                    auto& j1 = (el1.first)[1];
                    T& val1 = el1.second;
                    for (std::size_t j2 = 0; j2 < mat2.cols ; j2++)
                        
                        // Since mat2 is uncompressed, we use tha call operator in mat2 (find method).
                        result(i1,j2) += val1 * mat2(j1, j2);  
                }
            }

            if(mat1.compressed && !mat2.compressed){

                //mat1 is compressed, we need to save indices with the while method (same of compress()).
                std::size_t i1 = 0;
                for(std::size_t val_id1 = 0 ; val_id1 < mat1.values.size(); val_id1++){
                    while(mat1.row_indices[i1+1] < val_id1 + 1){
                        i1++; 
                    }
                    auto& j1 = mat1.col_indices[val_id1];
                    const T& val1 = mat1.values[val_id1];
                    for (std::size_t j2 = 0; j2 < mat2.cols ; j2++)
                        // Since mat2 is uncompressed, we use tha call operator in mat2 (find method).
                        result(i1,j2) += val1 * mat2(j1, j2);
                }
            }

            if(!mat1.compressed && mat2.compressed){
                
                //Same as above but the loops are inverted, to avoid initializing i2 multiple times.
                std::size_t i2 = 0;
                for(std::size_t val_id2 = 0 ; val_id2 < mat2.values.size(); val_id2++){
                    while(mat2.row_indices[i2+1] < val_id2 + 1){
                        i2++; 
                    }
                    const T& val2 = mat2.values[val_id2];
                    std::size_t j2 = mat2.col_indices[val_id2];
                    for (std::size_t i1 = 0; i1 < mat1.rows ; i1++){
                        const T& val1 = mat1(i1, i2);
                        result(i1,j2) += val1 * val2;
                    }
                }
            }

            // Only case in which it's not convenient to use the call operator. An if that checks 
            // if (j1==i2) is instead used.
            if(mat1.compressed && mat2.compressed){
                std::size_t i1 = 0;
                for(std::size_t val_id1 = 0 ; val_id1 < mat1.values.size(); val_id1++){
                    while(mat1.row_indices[i1+1] < val_id1 + 1){
                        i1++; 
                    }
                    std::size_t j1 = mat1.col_indices[val_id1];
                    const T& val1 = mat1.values[val_id1];
                    std::size_t i2 = 0;
                    for(std::size_t val_id2 = 0 ; val_id2 < mat2.values.size(); val_id2++){
                        while(mat2.row_indices[i2+1] < val_id2 + 1){
                            i2++; 
                        }
                        std::size_t j2 = mat2.col_indices[val_id2];
                        const T& val2 = mat2.values[val_id2];
                        if (j1 == i2){
                            result(i1,j2) += val1 * val2;
                        }
                    }
                }
            }

        }else if constexpr (Order1 == StorageOrder::COLUMN_MAJOR && Order2 == StorageOrder::COLUMN_MAJOR){
            
            // Same as above, but rows and cols indices are swapped in uncompressed format, as
            // well as the roles of col_indices and row_indices in compressed format.  
            // When both are compressed, the weighted sum of columns is implemented 
            // (like in matrix-vector product) .

            if(!mat1.compressed && !mat2.compressed){
                for(auto el1 : mat1.data){
                    auto& j1 = (el1.first)[0];
                    auto& i1 = (el1.first)[1];
                    T& val1 = el1.second;
                    for (std::size_t j2 = 0; j2 < mat2.cols ; j2++)
                        result(i1,j2) += val1 * mat2(j1, j2);
                }
            }

            if(mat1.compressed && !mat2.compressed){

                //Weighted sum of columns
                for (std::size_t j1 = 0; j1 < mat1.cols; ++j1) {
                    for (std::size_t i1 = mat1.col_indices[j1]; i1 < mat1.col_indices[j1 + 1]; ++i1) {
                        for(std::size_t j2 = 0 ; j2 < mat2.values.size() ; j2++){
                                result(mat1.row_indices[i1], j2) += mat1.values[i1] * mat2(j1,j2); 
                        }                                                                  
                    }
                }
            }

            if(!mat1.compressed && mat2.compressed){

                std::size_t j2 = 0;
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size() ; val_id2++){
                    while(mat2.col_indices[j2+1] < val_id2 + 1){
                            j2++; 
                        }
                    auto& i2 = mat2.row_indices[val_id2]; // ==j1
                    const T& val2 = mat2.values[val_id2];
                    for(std::size_t i1 = 0; i1 < mat1.rows; ++i1){
                        T val1 = mat1(i1, i2);
                        result(i1,j2) += val1 * val2;
                    }
                }
            }

            if(mat1.compressed && mat2.compressed){
                
                std::size_t j2 = 0;
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size(); val_id2++){
                    while(mat2.col_indices[j2+1] < val_id2 + 1){
                            j2++; 
                        }
                    auto& i2 = mat2.row_indices[val_id2];
                    const T& val2 = mat2.values[val_id2];
                    std::size_t j1 = 0;
                    for(std::size_t val_id1 = 0; val_id1 < mat1.values.size(); val_id1++){
                        while(mat1.col_indices[j1+1] < val_id1 + 1){
                            j1++; 
                        }
                        auto& i1 = mat1.row_indices[val_id1];
                        T val1 = mat1.values[val_id1];
                        if(j1 == i2){
                            result(i1,j2) += val1 * val2;
                        }
                        
                    }

                }
            }

        }else if constexpr (Order1 == StorageOrder::ROW_MAJOR && Order2 == StorageOrder::COLUMN_MAJOR){

            //All the below cases are repetitions of the above when needed.

            if(!mat1.compressed && !mat2.compressed){
                for(auto el1 : mat1.data){
                    auto& i1 = (el1.first)[0];
                    auto& j1 = (el1.first)[1];
                    T& val1 = el1.second;
                    for (std::size_t j2 = 0; j2 < mat2.cols ; j2++)
                        result(i1,j2) += val1 * mat2(j1, j2);
                }
            }

            if(mat1.compressed && !mat2.compressed){

                std::size_t i1 = 0;
                for(std::size_t val_id1 = 0; val_id1 < mat1.values.size(); ++val_id1) {
                    while(mat1.row_indices[i1+1] < val_id1 + 1){
                            i1++; 
                        }
                    auto& j1 = mat1.col_indices[val_id1];
                    const T& val1 = mat1.values[val_id1];
                    for (std::size_t j2 = 0; j2 < mat2.cols ; j2++)
                        result(i1,j2) += val1 * mat2(j1, j2);
                }
                    
            }

            if(!mat1.compressed && mat2.compressed){

                std::size_t j2 = 0;
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size() ; val_id2++){
                    while(mat2.col_indices[j2+1] < val_id2 + 1){
                            j2++; 
                        }
                    auto& i2 = mat2.row_indices[val_id2]; 
                    const T& val2 = mat2.values[val_id2];
                    for(std::size_t i1 = 0; i1 < mat1.rows; ++i1){
                        T val1 = mat1(i1, i2) ;
                        result(i1,j2) += val1 * val2;
                    }
                }
            }

            if(mat1.compressed && mat2.compressed){
                
                std::size_t j2 = 0;
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size(); val_id2++){
                    while(mat2.col_indices[j2+1] < val_id2 + 1){
                            j2++; 
                        }
                    auto& i2 = mat2.row_indices[val_id2];
                    const T& val2 = mat2.values[val_id2];
                    std::size_t i1 = 0;
                    for(std::size_t val_id1 = 0 ; val_id1 < mat1.values.size(); val_id1++){
                        while(mat1.row_indices[i1+1] < val_id1 + 1){
                            i1++; 
                        }
                        std::size_t j1 = mat1.col_indices[val_id1];
                        const T& val1 = mat1.values[val_id1];
                        if (j1 == i2)
                            result(i1,j2) += val1 * val2;
                    }

                }
            }

            

        }else if constexpr (Order1 == StorageOrder::COLUMN_MAJOR && Order2 == StorageOrder::ROW_MAJOR){

            if(!mat1.compressed && !mat2.compressed){
                for(auto el1 : mat1.data){
                    auto& j1 = (el1.first)[0];
                    auto& i1 = (el1.first)[1];
                    T& val1 = el1.second;
                    for (std::size_t j2 = 0; j2 < mat2.cols ; j2++)
                        result(i1,j2) += val1 * mat2(j1, j2);
                }
            }

            if(mat1.compressed && !mat2.compressed){

                for (std::size_t j1 = 0; j1 < mat1.cols; ++j1) {
                    for (std::size_t i1 = mat1.col_indices[j1]; i1 < mat1.col_indices[j1 + 1]; ++i1) {
                        for(std::size_t j2 = 0 ; j2 < mat2.values.size() ; j2++){
                                result(mat1.row_indices[i1], j2) += mat1.values[i1] * mat2(j1,j2); 
                        }                                                                  
                    }
                }
            }

            if(!mat1.compressed && mat2.compressed){
                
                std::size_t i2 = 0;
                for(std::size_t val_id2 = 0 ; val_id2 < mat2.values.size(); val_id2++){
                    while(mat2.row_indices[i2+1] < val_id2 + 1){
                        i2++; 
                    }
                    auto& j2 = mat2.col_indices[val_id2];
                    const T& val2 = mat2.values[val_id2];
                    for (std::size_t i1 = 0; i1 < mat1.rows ; i1++){
                        const T& val1 = mat1(i1, i2);
                        result(i1, j2) += val1 * val2;
                    }
                }
            } 

            if(mat1.compressed && mat2.compressed){
                
                std::size_t i2 = 0;
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size(); val_id2++){
                    while(mat2.row_indices[i2+1] < val_id2 + 1){
                        i2++; 
                    }
                    auto& j2 = mat2.col_indices[val_id2];
                    const T& val2 = mat2.values[val_id2];
                    std::size_t j1 = 0;
                    for(std::size_t val_id1 = 0; val_id1 < mat1.values.size(); val_id1++){
                        while(mat1.col_indices[j1+1] < val_id1 + 1){
                            j1++; 
                        }
                        auto& i1 = mat1.row_indices[val_id1];
                        T val1 = mat1.values[val_id1];
                        if(j1 == i2){
                            result(i1, j2) += val1 * val2;
                        }
                        
                    }

                }
            }

        }

        return result;    
    }

    template Matrix<double, StorageOrder::ROW_MAJOR> operator*<double, StorageOrder::ROW_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::ROW_MAJOR>& mat1, const Matrix<double, StorageOrder::ROW_MAJOR>& mat2);
    template Matrix<double, StorageOrder::COLUMN_MAJOR> operator*<double, StorageOrder::COLUMN_MAJOR, StorageOrder::COLUMN_MAJOR>(const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat1, const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat2);
    template Matrix<double, StorageOrder::ROW_MAJOR> operator*<double, StorageOrder::ROW_MAJOR, StorageOrder::COLUMN_MAJOR>(const Matrix<double, StorageOrder::ROW_MAJOR>& mat1, const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat2);
    template Matrix<double, StorageOrder::COLUMN_MAJOR> operator*<double, StorageOrder::COLUMN_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat1, const Matrix<double, StorageOrder::ROW_MAJOR>& mat2);
    
    template Matrix<std::complex<double>, StorageOrder::ROW_MAJOR> operator*<std::complex<double>, StorageOrder::ROW_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>& mat1, const Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>& mat2);
    template Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR> operator*<std::complex<double>, StorageOrder::COLUMN_MAJOR, StorageOrder::COLUMN_MAJOR>(const Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>& mat1, const Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>& mat2);
    template Matrix<std::complex<double>, StorageOrder::ROW_MAJOR> operator*<std::complex<double>, StorageOrder::ROW_MAJOR, StorageOrder::COLUMN_MAJOR>(const Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>& mat1, const Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>& mat2);
    template Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR> operator*<std::complex<double>, StorageOrder::COLUMN_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>& mat1, const Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>& mat2);

    
}