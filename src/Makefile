# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -I../include -O3 #-Werror
$(LINK.o)=CXX $(CXXFLAGS) $(CXXFLAGS) $(TARGET_ARCH)
# Source files
SRCS= MatMatProduct.cpp Product.cpp Norm.cpp Matrix.cpp 
LIB=libmymat.a
LDLIBS=-lmymat

#SRCS += main.cpp
#SRCS += mainComplex.cpp
#SRCS += mainMatMatProduct.cpp
# Phony targets
.PHONY: all clean

# Object files
OBJS= $(SRCS:.cpp=.o)
$(LIB): $(OBJS)
	ar rcs $@ $(OBJS)

# Executable name
EXEC=main#
EXEC_SRC=$(EXEC).cpp
EXEC_OBJ=$(EXEC).o
# Default target

all: $(EXEC)

comparison:
	$(MAKE) all EXEC=mainComparison

complex:
	$(MAKE) all EXEC=mainComplex

standard:
	$(MAKE) all EXEC=main

# Linking
$(EXEC): $(LIB) $(EXEC_OBJ) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(EXEC_OBJ) -L. $(LDLIBS)


# Clean
clean:
	rm -f $(OBJS) $(EXEC)

