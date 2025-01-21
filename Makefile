# Compiler
CXX = g++

# Compiler flags
INCLUDE = -I include/
CXX_FLAGS = -std=c++14

# Source file
SRCS = src/*.cpp

# Executable name
EXEC = build/h_p_res

# Build rule
$(EXEC): $(SRCS)
	$(CXX) $(GUROBI) $(INCLUDE) $(CXX_FLAGS) $(SRCS) -o $(EXEC)

# Clean rule
clean:
	rm -f $(EXEC)
