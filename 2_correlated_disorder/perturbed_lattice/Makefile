# Makefile for percolation simulation
# Optimized version without OpenMP

CC = gcc
CFLAGS = -Wall -Wextra -O3 -march=native -ffast-math -funroll-loops -ftree-vectorize
LDFLAGS = -lm

# Check for AVX2 support
AVX2_SUPPORT := $(shell echo | $(CC) -dM -E - -march=native 2>/dev/null | grep -c AVX2)
ifeq ($(AVX2_SUPPORT), 1)
    CFLAGS += -mavx2
endif

# Source files
SRCS = main.c lattice.c output.c percolation.c random_utils.c union_find.c

# Object files
OBJS = $(SRCS:.c=.o)

# Headers
HEADERS = config.h lattice.h output.h percolation.h random_utils.h union_find.h

# Target
TARGET = percolation

.PHONY: all clean debug profile test help

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Pattern rule for object files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Debug build
debug: CFLAGS = -Wall -Wextra -g -O0 -DDEBUG -fsanitize=address -fsanitize=undefined
debug: LDFLAGS += -fsanitize=address -fsanitize=undefined
debug: clean $(TARGET)

# Profile build (for gprof)
profile: CFLAGS += -pg
profile: LDFLAGS += -pg
profile: clean $(TARGET)

# Clean
clean:
	rm -f $(OBJS) $(TARGET) *.dat gmon.out

# Test with different weight functions
test: $(TARGET)
	@echo "Testing Power Law weight..."
	./$(TARGET) -L 50 -s 0.1 -w power -a 2.0 -T 5.0 -M 10 -D 10 -r 12345
	@echo "\nTesting Exponential weight..."
	./$(TARGET) -L 50 -s 0.1 -w exp -l 1.5 -T 5.0 -M 10 -D 10 -r 12345
	@echo "\nTesting Bernoulli weight..."
	./$(TARGET) -L 50 -s 0.1 -w bernoulli -l 1.0 -O 0.5 -T 5.0 -M 10 -D 10 -r 12345

# Benchmark
benchmark: $(TARGET)
	@echo "Benchmarking with L=100..."
	@echo "\nPower Law:"
	time ./$(TARGET) -L 100 -s 0.1 -w power -a 2.0 -T 10.0 -M 100 -D 100 -r 42
	@echo "\nExponential:"
	time ./$(TARGET) -L 100 -s 0.1 -w exp -l 1.5 -T 10.0 -M 100 -D 100 -r 42
	@echo "\nBernoulli:"
	time ./$(TARGET) -L 100 -s 0.1 -w bernoulli -l 1.0 -O 0.5 -T 10.0 -M 100 -D 100 -r 42

# Help
help:
	@echo "Available targets:"
	@echo "  all       - Build the percolation simulation"
	@echo "  debug     - Build with debug symbols and sanitizers"
	@echo "  profile   - Build with profiling support (gprof)"
	@echo "  test      - Run small tests with all weight functions"
	@echo "  benchmark - Run benchmarks with all weight functions"
	@echo "  clean     - Remove all built files"
	@echo "  help      - Show this help message"
	@echo ""
	@echo "Usage examples:"
	@echo "  Power law:    ./percolation -L 100 -s 0.1 -w power -a 2.0 -r 12345"
	@echo "  Exponential:  ./percolation -L 100 -s 0.1 -w exp -l 1.5 -r 12345"
	@echo "  Bernoulli:    ./percolation -L 100 -s 0.1 -w bernoulli -l 1.0 -O 0.5 -r 12345"
