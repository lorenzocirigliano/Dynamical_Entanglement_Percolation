CC = gcc
CFLAGS = -O3 -Wall -Wextra -std=c99 -march=native
LDFLAGS = -lm

# Source files
SRCS = main.c lattice.c percolation.c utils.c io.c
OBJS = $(SRCS:.c=.o)
TARGET = percolation_sim

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Rebuild everything
rebuild: clean all

# Run example simulation with alternating mode
run_example:
	./$(TARGET) alternating 50 100 0.5 0.0 1.0 0.0 1.0 0.1 12345

# Run quick test
test:
	./$(TARGET) alternating 10 10 0.5 0.5 0.5 0.5 0.5 0.1 42

.PHONY: all clean rebuild run_example test
