# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -I../include -O3 #-Werror

# Source files
SRCS := MatMatProduct.cpp Product.cpp Norm.cpp Matrix.cpp 
SRCS += mainComparison.cpp
#SRCS += main.cpp
#SRCS += mainComplex.cpp
#SRCS += mainMatMatProduct.cpp

# Object files
OBJS := $(SRCS:.cpp=.o)

# Executable name
EXEC := main 

# Default target
all: $(EXEC)

# Linking
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

# Compilation
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(EXEC)

# Phony targets
.PHONY: all clean
