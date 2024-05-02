# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -I../include #-Werror

# Source files
SRCS := Product.cpp Norm.cpp Matrix.cpp main.cpp

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
