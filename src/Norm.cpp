
#include "Matrix.hpp"
#include <cmath>
#include <algorithm>
#include <complex>

namespace algebra {    
    
    template<typename T, StorageOrder Order>
    template<NormType Norm>
    double Matrix<T, Order>::norm(){
        double ret = 0;
        if constexpr (Norm == NormType::One){
            if constexpr (Order == StorageOrder::ROW_MAJOR){
                std::vector<T> sums_along_columns(cols, 0);
                if(compressed){
                    for (std::size_t i = 0; i < rows; ++i) {
                        for (std::size_t j = row_indices[i]; j < row_indices[i + 1]; ++j) {
                            sums_along_columns[col_indices[j]] += std::abs(values[j]); 
                        }
                    }
                    return *std::max_element(sums_along_columns.begin(), 
                                            sums_along_columns.end(), 
                                            [](T a, T b){
                                                return std::abs(a) < std::abs(b);
                                            });
                    //return *ret;
                }else{
                    for (auto val : data){
                        auto& j = (val.first)[1];
                        sums_along_columns[j] += std::abs(val.second); 
                    }
                    return *std::max_element(sums_along_columns.begin(), 
                                            sums_along_columns.end(), 
                                            [](T a, T b){
                                                return std::abs(a) < std::abs(b);
                                            });
                    //return *ret;
                }
            }else{
                //std::vector<T> sums_along_columns(mat.cols, 0);
                T max_val = 0;
                if(compressed){
                    for (std::size_t j = 0; j < cols; ++j) {
                        T sum = 0;
                        for (std::size_t i = col_indices[j]; i < col_indices[j + 1]; ++i) {
                            //sums_along_columns[col_indices[j]] += std::abs(values[j]); 
                            sum += std::abs(values[i]);
                        }
                        max_val = std::max(max_val, sum);
                    }
                    return max_val;
                }else{
                    std::vector<T> sums_along_columns(cols, 0);
                    for (auto val : data){
                        auto& j = (val.first)[0];
                        sums_along_columns[j] += std::abs(val.second); 
                    }
                    return *std::max_element(sums_along_columns.begin(), 
                                            sums_along_columns.end(), 
                                            [](T a, T b){
                                                return std::abs(a) < std::abs(b);
                                            });
                    //return *ret;
                }
            }
        }else if constexpr(Norm == NormType::Infinity){
            if constexpr (Order == StorageOrder::ROW_MAJOR){
                if(compressed){
                    T max_val = 0;
                    for (std::size_t i = 0; i < rows; ++i) {
                        T sum = 0;
                        for (std::size_t j = row_indices[i]; j < row_indices[i + 1]; ++j) {
                            sum += std::abs(values[j]); 
                        }
                        max_val = std::max(max_val, sum);
                    }
                    return max_val;
                }else{
                    std::vector<T> sums_along_rows(rows, 0);
                    for (auto val : data){
                        auto& i = (val.first)[0];
                        sums_along_rows[i] += std::abs(val.second); 
                    }
                    return *std::max_element(sums_along_rows.begin(), 
                                            sums_along_rows.end(), 
                                            [](T a, T b){
                                                return std::abs(a) < std::abs(b);
                                            });
                    //return *ret;
                }
            }else{
                std::vector<T> sums_along_rows(rows, 0);
                if(compressed){
                    for (std::size_t j = 0; j < cols; ++j) {
                        for (std::size_t i = col_indices[j]; i < col_indices[j + 1]; ++i) {
                            sums_along_rows[row_indices[i]] += std::abs(values[i]); 
                        }
                    }
                    return *std::max_element(sums_along_rows.begin(), 
                                            sums_along_rows.end(), 
                                            [](T a, T b){
                                                return std::abs(a) < std::abs(b);
                                            });
                    //return *ret;
                }else{
                    std::vector<T> sums_along_rows(rows, 0);
                    for (auto val : data){
                        auto& i = (val.first)[1];
                        sums_along_rows[i] += std::abs(val.second); 
                    }
                    return *std::max_element(sums_along_rows.begin(), 
                                            sums_along_rows.end(), 
                                            [](T a, T b){
                                                return std::abs(a) < std::abs(b);
                                            });
                    //return *ret;
                }
            }
        }else if constexpr(Norm == NormType::Frobenius){
            if(compressed){
                for (auto val : values)
                    ret += val * val;
                return std::sqrt(ret);
            }else{
                for (auto val : data)
                    ret += val.second * val.second;
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
    
} // namespace algebra