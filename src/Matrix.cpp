
#include "Matrix.hpp"

namespace algebra {

    // Default constructor
    template<typename T, StorageOrder Order>
    Matrix<T, Order>::Matrix() : compressed(false), rows(0), cols(0) {}

    // Constructor with size
    template<typename T, StorageOrder Order>
    Matrix<T, Order>::Matrix(std::size_t size) : Matrix(size, size) {}

    // Constructor with rows and cols
    template<typename T, StorageOrder Order>
    Matrix<T, Order>::Matrix(std::size_t rows, std::size_t cols)
        : compressed(false), rows(rows), cols(cols) {}

    // Resize method
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::resize(std::size_t rows, std::size_t cols) {
        this->rows = rows;
        this->cols = cols;
        row_indices.clear();
        col_indices.clear();
        values.clear();
        data.clear();
        compressed = false;
    }

    // Compression method using CSR format
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::compress() {

        if (compressed) {
            return; // Matrix is already compressed
        }
        if constexpr (Order == StorageOrder::ROW_MAJOR){    

            // Fill CSR vectors from COOmap
            row_indices.resize(rows + 1); //first element is automatically set to 0.
            std::size_t count = 0;        //count of non-zero elements in current row 
            std::size_t row_idx = 0;      //last row index

            for (auto& entry : data) {

                //extracting current element index and value
                auto& current_row = (entry.first)[0];
                auto& current_col = (entry.first)[1];
                T value = entry.second;

                //updating the current row index
                while(row_idx < current_row){
                    /*  
                    When we enter a new row, the new row index is saved as previous_row_index + count, 
                    then count is set to 0, in this way we solve the case of empty rows
                    */
                    row_idx++; 
                    row_indices[row_idx] = row_indices[row_idx - 1] + count;
                    count = 0;
                }
                col_indices.push_back(current_col);
                values.push_back(value);
                ++count;  
            } 
            row_indices[row_idx + 1] = row_indices[row_idx] + count;

        }else{

            // Fill CSC vectors from COOmap
            // Identical to the above case but with roles of row_indices and column_indeces swapped

            col_indices.resize(cols + 1); //first element is automatically set to 0.
            std::size_t count = 0;        //count of non-zero elements in current column 
            std::size_t col_idx = 0;      //last column index

            for (const auto& entry : data) {
                auto& current_col = (entry.first)[0];
                auto& current_row = (entry.first)[1];
                T value = entry.second;
                while(col_idx < current_col){
                    col_idx++; 
                    col_indices[col_idx] = col_indices[col_idx - 1] + count;
                    count = 0;
                }
                row_indices.push_back(current_row);
                values.push_back(value);
                ++count;  
            } 
            col_indices[col_idx + 1] = col_indices[col_idx] + count;
        }

        // Clear COOmap
        data.clear();

        // Set compressed flag
        compressed = true;
    }

    // Uncompression method
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::uncompress() {

        if (!compressed) {
            return; // Matrix is already uncompressed
        }
        
        //Selecting elements one row at a time
        if constexpr (Order==StorageOrder::ROW_MAJOR){

            //ROW_MAJOR

            for(std::size_t i = 0; i < rows; ++i){

                //Number of elements in row i = row_indices[i+1] - row_indices[i] 
                for (std::size_t j = row_indices[i]; j < row_indices[i+1]; ++j){
                    data[{i, col_indices[j]}] = values[j];
                };
            }
        }else{
            
            //COLUMN_MAJOR
                 
            for(std::size_t j = 0; j < cols; ++j){

                //Number of elements in col i = col_indices[i+1] - col_indices[i]
                for (std::size_t i = col_indices[j]; i < col_indices[j+1]; ++i){

                    //rows and columns are swapped in the key array
                    data[{j, row_indices[i]}] = values[i];
                };
            }
            
        }

        // Clear CSR vectors
        row_indices.clear();
        col_indices.clear();
        values.clear();

        // Set compressed flag
        compressed = false;
    }

    // Check if matrix is compressed
    template<typename T, StorageOrder Order>
    bool Matrix<T, Order>::isCompressed() const {
        return compressed;
    }

    // Const version of operator()
    template<typename T, StorageOrder Order>
    T Matrix<T, Order>::operator()(std::size_t i, std::size_t j) const {
        if (i >= rows || j >= cols) {
            throw std::out_of_range("Matrix indices out of bounds");
        }
        if constexpr (Order==StorageOrder::ROW_MAJOR){
            if (compressed) {

                // Search for element in CSR format 
                // n_elements_in_row = row_indices[i+1] - row_indices[i];
                for (std::size_t k = row_indices[i]; k < row_indices[i+1]; ++k){
                    if(col_indices[k] == j){
                        return values[k];
                        }
                    }
                return 0; // Element not found, default value is 0

            } else {

                // Search for element in COOmap
                auto it = data.find({i, j});
                if (it == data.end()) {
                    return 0; // Element not found, default value is 0
                }
                return it->second;
            }
        }else{
            if (compressed) {

                // Search for element in CSC format

                for (std::size_t k = col_indices[j]; k < col_indices[j+1]; ++k){
                    if(row_indices[k] == i){
                        return values[k];
                        }
                    }
                //std::cout<<"Executing compressed () with ColumnMajor ordering"<<std::endl;
                return 0; // Element not found, assuming default value is 0

            } else {

                // Search for element in COOmap
                auto it = data.find({j, i});
                if (it == data.end()) {
                    return 0; // Element not found, default value is 0
                }
                return it->second;
            }
        }
    }

    // Non-const version of operator()
    template<typename T, StorageOrder Order>
    T& Matrix<T, Order>::operator()(std::size_t i, std::size_t j) {
        if (i >= rows || j >= cols) {
            throw std::out_of_range("Matrix indices out of bounds");
        }
        if constexpr (Order == StorageOrder::ROW_MAJOR){
            if (compressed) {

                // Search for element in CSR format 
                // n_elements_in_row = row_indices[i+1] - row_indices[i];
                for (std::size_t k = row_indices[i]; k < row_indices[i+1]; ++k){
                    if(col_indices[k] == j)
                        return values[k];
                }

               throw std::runtime_error("Element not found, switch to uncompressed format to add it."); 

            } else {
                // Search for element in COOmap
                return data[{i, j}];
            }
        }else{
            if (compressed) {

                // Search for element in CSR format (my version)
                // n_elements_in_row = row_indices[i+1] - row_indices[i];
                for (std::size_t k = col_indices[j]; i < col_indices[j+1]; ++k){
                    if(row_indices[k] == i)
                        return values[k];
                }

                throw std::runtime_error("Element not found, switch to uncompressed format to add it.");

            } else {
                // Search for element in COOmap
                return data[{j, i}];
            }
        }
    }

    
    
    // Explicit instantiations
    template class Matrix<int, StorageOrder::ROW_MAJOR>;
    template class Matrix<double, StorageOrder::ROW_MAJOR>;
    template class Matrix<float, StorageOrder::ROW_MAJOR>;
    template class Matrix<std::complex<double>, StorageOrder::ROW_MAJOR>;

    template class Matrix<int, StorageOrder::COLUMN_MAJOR>;
    template class Matrix<double, StorageOrder::COLUMN_MAJOR>;
    template class Matrix<float, StorageOrder::COLUMN_MAJOR>;
    template class Matrix<std::complex<double>, StorageOrder::COLUMN_MAJOR>;

} // namespace algebra