#pragma once

#include <stdlib.h>
#include <regex.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include "structure.h"

#define TABLE_SIZE 100000

typedef struct entry {
    int* key;
    int key_size;
    double value;
    struct entry* next;
} entry;

typedef struct {
    entry* buckets[TABLE_SIZE];
} hash_table;

int hash_function(int* key, int key_size);
hash_table* create_hash_table();
int arrays_equal(int* arr1, int* arr2, int size);
void insert(hash_table* table, int* key, int key_size, double value);
double search(hash_table* table, int* key, int key_size);
void delete_hash_table(hash_table* table);
int _exists_in_subarray(int* arr, int start, int end, int value);
int _min_value_index(double* values, size_t size);
int compare(const void* a, const void* b);
void _reverse_subarray(int* arr, int start, int end);
int compare_fitness(const void* a, const void* b);
void write_array_to_csv(const char* filename, double** arr, size_t size_x, size_t size_y);
void write_string_array_to_csv(const char* filename, char*** arr, size_t size_x, size_t size_y);
void write_string_array_to_csv(const char* filename, char*** arr, size_t size_x, size_t size_y);
void save_result(invidual*** arr, size_t size_x, size_t size_y, const char* filename);