
#include "Matrix.hpp"

namespace algebra{

    //Product with general matrices
    template<typename T, StorageOrder Order1, StorageOrder Order2, StorageOrder OrderOut>
    Matrix<T, OrderOut> operator*(const Matrix<T, Order1>& mat1, const Matrix<T, Order2>& mat2) {

        if (mat2.rows != mat1.cols) 
            throw std::out_of_range("Incompatible sizes");    

        Matrix<T, OrderOut> result(mat1.rows, mat2.cols);

        if constexpr (Order1 == StorageOrder::ROW_MAJOR && Order2 == StorageOrder::ROW_MAJOR){

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
                for(std::size_t val_id1 = 0 ; val_id1 < mat1.values.size(); val_id1++){
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
                
                //alternative impl
                std::size_t i2 = 0;
                for(std::size_t val_id2 = 0 ; val_id2 < mat2.values.size(); val_id2++){
                    while(mat2.row_indices[i2+1] < val_id2 + 1){
                        i2++; 
                    }
                    const T& val2 = mat2.values[val_id2];
                    std::size_t j2 = mat2.col_indices[val_id2];
                    for (std::size_t i1 = 0; i1 < mat1.rows ; i1++){
                        //i2 == j1
                        const T& val1 = mat1(i1, i2);// ==0 the call operator gives zero if element is not present
                        result(i1,j2) += val1 * val2;
                    }
                }
            }

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

            //throw std::runtime_error("Not implemented yet");

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
                                result(mat1.row_indices[i1], j2) += mat1.values[i1] * mat2(j1,j2); //is there a better way? 
                        }                                                                  //Should work because mat2 is uncompressed
                    }
                }
            }

            if(!mat1.compressed && mat2.compressed){

                //alternative impl
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

            //throw std::runtime_error("Not implemented yet");

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

                //alternative impl
                std::size_t j2 = 0;
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size() ; val_id2++){
                    while(mat2.col_indices[j2+1] < val_id2 + 1){
                            j2++; 
                        }
                    auto& i2 = mat2.row_indices[val_id2]; // ==j1
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

                //Weighted sum of columns
                for (std::size_t j1 = 0; j1 < mat1.cols; ++j1) {
                    for (std::size_t i1 = mat1.col_indices[j1]; i1 < mat1.col_indices[j1 + 1]; ++i1) {
                        for(std::size_t j2 = 0 ; j2 < mat2.values.size() ; j2++){
                                result(mat1.row_indices[i1], j2) += mat1.values[i1] * mat2(j1,j2); //is there a better way? 
                        }                                                                  //Should work because mat2 is uncompressed
                    }
                }
            }

            if(!mat1.compressed && mat2.compressed){
                
                //alternative impl
                std::size_t i2 = 0;
                for(std::size_t val_id2 = 0 ; val_id2 < mat2.values.size(); val_id2++){
                    while(mat2.row_indices[i2+1] < val_id2 + 1){
                        i2++; 
                    }
                    auto& j2 = mat2.col_indices[val_id2];
                    const T& val2 = mat2.values[val_id2];
                    for (std::size_t i1 = 0; i1 < mat1.rows ; i1++){
                        //i2 == j1
                        const T& val1 = mat1(i1, i2);// ==0 the call operator gives zero if element is not present
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

    template Matrix<double, StorageOrder::ROW_MAJOR> operator*<double, StorageOrder::ROW_MAJOR, StorageOrder::ROW_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::ROW_MAJOR>& mat1, const Matrix<double, StorageOrder::ROW_MAJOR>& mat2);
    template Matrix<double, StorageOrder::ROW_MAJOR> operator*<double, StorageOrder::COLUMN_MAJOR, StorageOrder::COLUMN_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat1, const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat2);
    template Matrix<double, StorageOrder::ROW_MAJOR> operator*<double, StorageOrder::ROW_MAJOR, StorageOrder::COLUMN_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::ROW_MAJOR>& mat1, const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat2);
    template Matrix<double, StorageOrder::ROW_MAJOR> operator*<double, StorageOrder::COLUMN_MAJOR, StorageOrder::ROW_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat1, const Matrix<double, StorageOrder::ROW_MAJOR>& mat2);

}