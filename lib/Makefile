CC = gcc
CFLAGS = -std=c99 -O3 -fPIC
OBJ = main.o algorithms.o utils.o structure.o
TARGET = metaheuristics.so

$(TARGET): $(OBJ)
	$(CC) -shared $(OBJ) -o $(TARGET) -lm

algorithms.o: algorithms.c algorithms.h utils.h structure.h
	$(CC) $(CFLAGS) -c algorithms.c

utils.o: utils.c utils.h algorithms.h structure.h
	$(CC) $(CFLAGS) -c utils.c

structure.o: structure.c structure.h algorithms.h utils.h
	$(CC) $(CFLAGS) -c structure.c

main.o: main.c algorithms.h structure.h utils.h
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -f $(OBJ)