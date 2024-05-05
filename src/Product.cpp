
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
    template<typename T, StorageOrder Order1, StorageOrder Order2>
    std::vector<T> operator*(const Matrix<T, Order1>& mat1, const Matrix<T, Order2>& mat2) {

        if (mat2.cols > 1) 
            throw std::out_of_range("Only 1-column matrices are accepted");

        if (mat2.rows != mat1.cols) 
            throw std::out_of_range("Incompatible sizes");    

        std::vector<T> result(mat1.rows, 0);

        if constexpr (Order1 == StorageOrder::ROW_MAJOR && Order2 == StorageOrder::ROW_MAJOR){

            if(!mat1.compressed && !mat2.compressed){
                for(auto el1 : mat1.data){
                    auto& i1 = (el1.first)[0];
                    auto& j1 = (el1.first)[1];
                    T& val1 = el1.second;
                    result[i1] += val1 * mat2(j1, 0);
                }
            }

            if(mat1.compressed && !mat2.compressed){
                std::size_t i1 = 0;
                for(std::size_t val_id1 = 0 ; val_id1 < mat1.values.size(); val_id1++){
                    while(mat1.row_indices[i1+1] < val_id1 + 1){
                        i1++; 
                    }
                    std::size_t j1 = mat1.col_indices[val_id1];
                    const T& val1 = mat1.values[val_id1];
                    result[i1] += val1 * mat2(j1, 0);
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
                    for (std::size_t i1 = 0; i1 < mat1.rows ; i1++){
                        //i2 == j1
                        const T& val1 = mat1(i1, i2);// ==0 the call operator gives zero if element is not present
                        result[i1] += val1 * val2;
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
                            result[i1] += val1 * val2;
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

                    result[i1] += val1 * mat2(j1, 0);
                }
            }

            if(mat1.compressed && !mat2.compressed){

                //Weighted sum of columns
                for (std::size_t j = 0; j < mat1.cols; ++j) {
                    for (std::size_t i = mat1.col_indices[j]; i < mat1.col_indices[j + 1]; ++i) {
                        result[mat1.row_indices[i]] += mat1.values[i] * mat2(j,0); //is there a better way? 
                                                                                   //Should work because mat2 is uncompressed
                    }
                }
            }

            if(!mat1.compressed && mat2.compressed){

                //alternative impl
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size() ; val_id2++){
                    auto& i2 = mat2.row_indices[val_id2]; // ==j1
                    const T& val2 = mat2.values[val_id2];
                    for(std::size_t i1 = 0; i1 < mat1.rows; ++i1){
                        T val1 = mat1(i1, i2) ;
                        result[i1] += val1 * val2;
                    }
                }
            }

            if(mat1.compressed && mat2.compressed){

                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size(); val_id2++){
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
                            result[i1] += val1 * val2;
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

                    result[i1] += val1 * mat2(j1, 0);
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
                    result[i1] += mat1.values[val_id1] * mat2(j1,0);
                }
                    
            }

            if(!mat1.compressed && mat2.compressed){

                //alternative impl
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size() ; val_id2++){
                    auto& i2 = mat2.row_indices[val_id2]; // ==j1
                    const T& val2 = mat2.values[val_id2];
                    for(std::size_t i1 = 0; i1 < mat1.rows; ++i1){
                        T val1 = mat1(i1, i2) ;
                        result[i1] += val1 * val2;
                    }
                }
            }

            if(mat1.compressed && mat2.compressed){

                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size(); val_id2++){
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
                            result[i1] += val1 * val2;
                    }

                }
            }

            

        }else if constexpr (Order1 == StorageOrder::COLUMN_MAJOR && Order2 == StorageOrder::ROW_MAJOR){

            if(!mat1.compressed && !mat2.compressed){
                for(auto el1 : mat1.data){
                    auto& j1 = (el1.first)[0];
                    auto& i1 = (el1.first)[1];
                    T& val1 = el1.second;

                    result[i1] += val1 * mat2(j1, 0);
                }
            }

            if(mat1.compressed && !mat2.compressed){

                //Weighted sum of columns
                for (std::size_t j = 0; j < mat1.cols; ++j) {
                    for (std::size_t i = mat1.col_indices[j]; i < mat1.col_indices[j + 1]; ++i) {
                        result[mat1.row_indices[i]] += mat1.values[i] * mat2(j,0); //is there a better way? 
                                                                                   //Should work because mat2 is uncompressed
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
                    const T& val2 = mat2.values[val_id2];
                    for (std::size_t i1 = 0; i1 < mat1.rows ; i1++){
                        //i2 == j1
                        const T& val1 = mat1(i1, i2);// ==0 the call operator gives zero if element is not present
                        result[i1] += val1 * val2;
                    }
                }
            } 

            if(mat1.compressed && mat2.compressed){
                
                std::size_t i2 = 0;
                for(std::size_t val_id2 = 0; val_id2 < mat2.values.size(); val_id2++){
                    while(mat2.row_indices[i2+1] < val_id2 + 1){
                        i2++; 
                    }
                    const T& val2 = mat2.values[val_id2];
                    std::size_t j1 = 0;
                    for(std::size_t val_id1 = 0; val_id1 < mat1.values.size(); val_id1++){
                        while(mat1.col_indices[j1+1] < val_id1 + 1){
                            j1++; 
                        }
                        auto& i1 = mat1.row_indices[val_id1];
                        T val1 = mat1.values[val_id1];
                        if(j1 == i2){
                            result[i1] += val1 * val2;
                        }
                        
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

    template std::vector<double> operator*<double, StorageOrder::ROW_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::ROW_MAJOR>& mat1, const Matrix<double, StorageOrder::ROW_MAJOR>& mat2);
    template std::vector<double> operator*<double, StorageOrder::COLUMN_MAJOR, StorageOrder::COLUMN_MAJOR>(const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat1, const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat2);
    template std::vector<double> operator*<double, StorageOrder::ROW_MAJOR, StorageOrder::COLUMN_MAJOR>(const Matrix<double, StorageOrder::ROW_MAJOR>& mat1, const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat2);
    template std::vector<double> operator*<double, StorageOrder::COLUMN_MAJOR, StorageOrder::ROW_MAJOR>(const Matrix<double, StorageOrder::COLUMN_MAJOR>& mat1, const Matrix<double, StorageOrder::ROW_MAJOR>& mat2);

    //All of this must be done with std::complex<double> instead of double
}