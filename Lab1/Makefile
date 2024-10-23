# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra -std=c99 -O3

# Target exec
TARGET = main

# Source files
SRC = main.c

# Object files
OBJ = $(SRC:.c=.o)

# Defualt target
all: $(TARGET)

# Rule to build target exec
$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

# Rule to compile source files into object files
%.o: %.c
	$(CC)  $(CFLAGS) -c $< -o $@

# Rule to clean
clean:
	rm -f $(OBJ) $(TARGET)

# Rule to run exec
run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run