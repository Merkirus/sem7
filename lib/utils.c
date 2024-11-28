#include "utils.h"
#include "structure.h"
#include <math.h>

int hash_function(int* key, int key_size) {
    int hash = 0;
    for (int i = 0; i < key_size; i++)
        hash = (hash * 31 + key[i]) % TABLE_SIZE;
    return hash;
}

hash_table* create_hash_table() {
    hash_table* table = (hash_table*)malloc(sizeof(hash_table));
    for (int i = 0; i < TABLE_SIZE; i++)
        table->buckets[i] = NULL;
    return table;
}

int arrays_equal(int* arr1, int* arr2, int size) {
    for (int i = 0; i < size; i++) {
        if (arr1[i] != arr2[i])
            return 0;
    }
    return 1;
}

void insert(hash_table* table, int* key, int key_size, double value) {
    int index = hash_function(key, key_size);
    entry* current = table->buckets[index];

    while (current) {
        if (current->key_size == key_size && arrays_equal(current->key, key, key_size)) {
            current->value = value;
            return;
        }
        current = current->next;
    }

    entry* new_entry = (entry*)malloc(sizeof(entry));
    new_entry->key = (int*)malloc(key_size * sizeof(int));
    memcpy(new_entry->key, key, key_size * sizeof(int));
    new_entry->key_size = key_size;
    new_entry->value = value;
    new_entry->next = table->buckets[index];
    table->buckets[index] = new_entry;
}

double search(hash_table* table, int* key, int key_size) {
    int index = hash_function(key, key_size);
    entry* current = table->buckets[index];

    while (current) {
        if (current->key_size == key_size && arrays_equal(current->key, key, key_size))
            return current->value;
        current = current->next;
    }
    return NAN;
}

void delete_hash_table(hash_table* table) {
    for (int i = 0; i < TABLE_SIZE; i++) {
        entry* current = table->buckets[i];
        while (current) {
            entry* temp = current;
            current = current->next;
            free(temp->key);
            free(temp);
        }
    }
    free(table);
}

int _exists_in_subarray(int* arr, int start, int end, int value) {
    for (int i = start; i <= end; i++) {
        if (arr[i] == value)
            return 1;
        return 0;
    }
}

int _min_value_index(double* values, size_t size) {
    int min = INT_MAX;
    int min_index = 0;

    for (int i = 0; i < size; i++) {
        int value = values[i];

        if (value != 0 && value < min) {
            min = value;
            min_index = i;
        }
    }

    return min_index;
}

int compare(const void* a, const void* b) {
    return (*(double*)a - *(double*)b);
}

void _reverse_subarray(int* arr, int start, int end) {
    while (start < end) {
        int temp = arr[start];
        arr[start] = arr[end];
        arr[end] = temp;
        start++;
        end--;
    }
}

int compare_fitness(const void* a, const void* b) {
    invidual* inv1 = *(invidual**)a;
    invidual* inv2 = *(invidual**)b;

    if (inv1->fitness < inv2->fitness) return -1;
    if (inv1->fitness > inv2->fitness) return 1;
    return 0;
}

void write_array_to_csv(const char* filename, double** arr, size_t size_x, size_t size_y) {
    FILE *file = fopen(filename, "w"); // Open the file for writing
    if (file == NULL) {
        printf("Error opening file %s for writing\n", filename);
        return;
    }

    // Loop through each row
    for (int i = 0; i < size_x; i++) {
        // Loop through each column in the row
        for (int j = 0; j < size_y; j++) {
            fprintf(file, "%f", arr[i][j]);
            if (j < size_y - 1) {
                fprintf(file, ","); // Add a comma after each value except the last one
            }
        }
        fprintf(file, "\n"); // Add a newline after each row
    }

    // Close the file
    fclose(file);
    printf("2D array written to %s successfully.\n", filename);
}


void write_string_array_to_csv(const char* filename, char*** arr, size_t size_x, size_t size_y) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file %s for writing\n", filename);
        return;
    }


    for (int i = 0; i < size_x; i++) {
        // Loop through each column in the row
        for (int j = 0; j < size_y; j++) {
            fprintf(file, "%s", arr[i][j]);
            if (j < size_y - 1) {
                fprintf(file, ","); // Add a comma after each value except the last one
            }
        }
        fprintf(file, "\n"); // Add a newline after each row
    }

    fclose(file);
    printf("2D array written to %s successfully.\n", filename);
}

void save_result(invidual*** arr, size_t x, size_t y, const char* filename) {
    double** csv_arr = (double**)malloc(x * sizeof(double*));
    for (int i = 0; i < x; i++)
        csv_arr[i] = (double*)malloc(y * sizeof(double));

    char*** csv_str_arr = (char***)malloc(x * sizeof(char**));
    for (int i = 0; i < x; i++)
        csv_str_arr[i] = (char**)malloc(y * sizeof(char*));
    
    size_t gene_size = (arr[0][0]->gene.size * 6) + 3; // and terminator

    for (int i = 0; i < x; i++)
        for (int j = 0; j < y; j++)
            csv_str_arr[i][j] = (char*)malloc(gene_size * sizeof(char));
    
    for (int i = 0; i < x; i++)
        for (int j = 0; j < y; j++)
            csv_arr[i][j] = arr[i][j]->fitness;
    
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            csv_str_arr[i][j][0] = '\0';
            strcat(csv_str_arr[i][j], "(");

            for (int k = 0; k < arr[i][j]->gene.size; k++) {
                char temp[10];
                snprintf(temp, sizeof(temp), "%d", arr[i][j]->gene.entries[k]);
                strncat(csv_str_arr[i][j], temp, gene_size - strlen(csv_str_arr[i][j]) - 1);
                if (k < arr[i][j]->gene.size - 1)
                    strncat(csv_str_arr[i][j], ";", gene_size - strlen(csv_str_arr[i][j]) - 1);
            }

            strcat(csv_str_arr[i][j], ")");
        }
    }
    
    size_t str_filename_size = strlen(filename);
    str_filename_size += 5;
    char* str_filename[str_filename_size + 1];
    snprintf(str_filename, sizeof(str_filename), "str_%s", filename);
    write_array_to_csv(filename, csv_arr, x, y);
    write_string_array_to_csv(str_filename, csv_str_arr, x, y);
    for (int i = 0; i < x; i++)
        free(csv_arr[i]);
    free(csv_arr);
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++)
            free(csv_str_arr[i][j]);
        free(csv_str_arr[i]);
    }
    free(csv_str_arr);
}