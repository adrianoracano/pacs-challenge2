
#include "Matrix.hpp"
#include <cmath>
#include <algorithm>
#include <type_traits>

namespace algebra {  

    template<typename T>
    T max_element(std::vector<T>& vec){
        return *std::max_element(vec.begin(), vec.end(), 
                                            [](T a, T b){
                                                return std::abs(a) < std::abs(b);
                                            });
    }  
    
    template<typename T, StorageOrder Order>
    template<NormType Norm>
    double Matrix<T, Order>::norm(){
        double ret = 0;
        if constexpr (Norm == NormType::One){
            if constexpr (Order == StorageOrder::ROW_MAJOR){
                std::vector<double> sums_along_columns(cols, 0);
                if(compressed){
                    for (std::size_t i = 0; i < rows; ++i) {
                        for (std::size_t j = row_indices[i]; j < row_indices[i + 1]; ++j) {
                            sums_along_columns[col_indices[j]] += std::abs(values[j]); 
                        }
                    }
                    return max_element(sums_along_columns);
                }else{
                    for (auto val : data){
                        auto& j = (val.first)[1];
                        sums_along_columns[j] += std::abs(val.second); 
                    }
                    return max_element(sums_along_columns);
                    //return *ret;
                }
            }else{
                //std::vector<double> sums_along_columns(mat.cols, 0);
                double max_val = 0;
                if(compressed){
                    for (std::size_t j = 0; j < cols; ++j) {
                        double sum = 0;
                        for (std::size_t i = col_indices[j]; i < col_indices[j + 1]; ++i) {
                            //sums_along_columns[col_indices[j]] += std::abs(values[j]); 
                            sum += std::abs(values[i]);
                        }
                        max_val = std::max(std::abs(max_val), std::abs(sum));
                    }
                    return max_val;
                }else{
                    std::vector<double> sums_along_columns(cols, 0);
                    for (auto val : data){
                        auto& j = (val.first)[0];
                        sums_along_columns[j] += std::abs(val.second); 
                    }
                    return max_element(sums_along_columns);
                    //return *ret;
                }
            }
        }else if constexpr(Norm == NormType::Infinity){
            if constexpr (Order == StorageOrder::ROW_MAJOR){
                if(compressed){
                    double max_val = 0;
                    for (std::size_t i = 0; i < rows; ++i) {
                        double sum = 0;
                        for (std::size_t j = row_indices[i]; j < row_indices[i + 1]; ++j) {
                            sum += std::abs(values[j]); 
                        }
                        max_val = std::max(max_val, sum);
                    }
                    return max_val;
                }else{
                    std::vector<double> sums_along_rows(rows, 0);
                    for (auto val : data){
                        auto& i = (val.first)[0];
                        sums_along_rows[i] += std::abs(val.second); 
                    }
                    return max_element(sums_along_rows);
                    //return *ret;
                }
            }else{
                std::vector<double> sums_along_rows(rows, 0);
                if(compressed){
                    for (std::size_t j = 0; j < cols; ++j) {
                        for (std::size_t i = col_indices[j]; i < col_indices[j + 1]; ++i) {
                            sums_along_rows[row_indices[i]] += std::abs(values[i]); 
                        }
                    }
                    return max_element(sums_along_rows);
                    //return *ret;
                }else{
                    std::vector<double> sums_along_rows(rows, 0);
                    for (auto val : data){
                        auto& i = (val.first)[1];
                        sums_along_rows[i] += std::abs(val.second); 
                    }
                    return max_element(sums_along_rows);
                    //return *ret;
                }
            }
        }else if constexpr(Norm == NormType::Frobenius){
            if(compressed){
                for (auto val : values)
                    ret += std::abs(val) * std::abs(val);
                return std::sqrt(ret);
            }else{
                for (auto val : data)
                    ret += std::abs(val.second) * std::abs(val.second);
                return std::sqrt(ret);
            }
        }
        
    };
    
    //Explicit instantiation for norm function
    template double Matrix<double, StorageOrder::ROW_MAJOR>::norm<NormType::One>() ;
    template double Matrix<double, StorageOrder::ROW_MAJOR>::norm<NormType::Frobenius>() ;
    template double Matrix<double, StorageOrder::ROW_MAJOR>::norm<NormType::Infinity>() ;
    template double Matrix<double, StorageOrder::COLUMN_MAJOR>::norm<NormType::One>() ;
    template double Matrix<double, StorageOrder::COLUMN_MAJOR>::norm<NormType::Frobenius>() ;
    template double Matrix<double, StorageOrder::COLUMN_MAJOR>::norm<NormType::Infinity>() ;
    template double Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>::norm<NormType::One>() ;
    template double Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>::norm<NormType::Frobenius>() ;
    template double Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>::norm<NormType::Infinity>() ;
    template double Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>::norm<NormType::One>() ;
    template double Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>::norm<NormType::Frobenius>() ;
    template double Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>::norm<NormType::Infinity>() ;
    
} // namespace algebra