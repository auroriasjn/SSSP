# -------------------------
# Compiler
# -------------------------
CXX := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -pedantic -I./parlaylib/include -pthread
CXXFLAGS += -DNDEBUG
# CXXFLAGS += -fopenmp
# LDFLAGS  := -fopenmp
LDFLAGS := -pthread

# -------------------------
# Directories
# -------------------------
BUILD_DIR := build
TARGET := $(BUILD_DIR)/sssp_solver

# -------------------------
# Source files
# -------------------------
SRCS := \
	sssp.cpp \
	sssp/dijkstra.cpp \
	sssp/parallel_dijkstra.cpp \
	sssp/bundle_dijkstra.cpp \
	sssp/parallel_bundle_dijkstra.cpp \
	sssp/rho_stepping.cpp \
	sssp/bellman_ford.cpp \
	benchmark/source_selector.cpp \
	benchmark/benchmark_solver.cpp \
	parser/parser.cpp

# Object files (mirror directory structure inside build/)
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)

# -------------------------
# Default target
# -------------------------
all: $(TARGET)

# -------------------------
# Link
# -------------------------
$(TARGET): $(OBJS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# -------------------------
# Compile
# -------------------------
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# -------------------------
# Clean
# -------------------------
clean:
	rm -rf $(BUILD_DIR)

rebuild: clean all

.PHONY: all clean rebuild