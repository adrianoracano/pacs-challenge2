# PACS Challenge (2):
### Brief explanation:
The above code provides the implementation of a sparse matrix class template . 
The main methods are: 
- compress(), uncompress(): shift from COOmap to compressed format;
- operator() : call operator to access or add elements;
- norm() : template function for the norm (can be One, Infinity or Frobenius);
- friend operator*() : matrix vector product;
- friend operator*() : matrix matrix product.

All the methods are available for the types `double`, `std::complex<double>` and both in row major and column major orderings.

### How to run the code:
The repository is divided in:
- include : (header files)
- src : (source files)
To run the program it is necessary to type `make` and `./main`. It is possible to execute other versions of `main` (`mainComplex.cpp`, `mainMatMatProduct.cpp`) by selecting them in the makefile.
