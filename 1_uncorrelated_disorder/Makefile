CC = gcc
CFLAGS = -Wall -O3 -std=c99
LDFLAGS = -lm

# Source files
SRCS = main.c lattice.c union_find.c percolation.c random_utils.c output.c
OBJS = $(SRCS:.c=.o)

# Executable name
TARGET = square_uncorrelated

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
