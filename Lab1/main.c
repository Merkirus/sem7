#pragma once

#include <stdlib.h> 
#include <stddef.h>
#include <math.h>
#include <regex.h>
#include <string.h>

#define RANDOM(min, max) (rand() % ((max) - (min) - (min) + 1) + (min))
#define BLANK -100

double _distance(double x1, double y1, double x2, double y2);
char* _get_file_data(char* line);
int _exists_in_subarray(int* arr, int start, int end, int value);

typedef struct {
    size_t size;
    int* entries;
} vec2d;

typedef struct {
    int size_x;
    int size_y;
    int** entries;
} vec3d;

typedef struct {
    vec2d gene;
} invidual;

typedef struct {
    int dim;
    int cap;
    double** entries;
    int* demand;
} problem;

invidual* create_invidual(size_t n);
void delete_invidual(invidual* inv);

problem* create_problem(char* file_path);

invidual* crossover_ox(invidual* inv1, invidual* inv2);
invidual* crossover_pmx(invidual* inv1, invidual* inv2);
invidual* crossover_cx(invidual* inv1, invidual* inv2);

double alg_greedy(problem* prob, invidual* inv);
void mutation_swap(invidual* inv, double prob);

double _distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
}

char* _get_file_data(char* line) {
    regex_t regex;
    regmatch_t matches[2];
    char* pattern = "\\w+\\s*:\\s*(.*)";
    char * result = NULL;
    
    if (regcomp(&regex, pattern, REG_EXTENDED) != 0)
        return NULL;

    if (regexec(&regex, line, 2, matches, 0) == 0) {
        if (matches[1].rm_so != -1) {
            int length = matches[1].rm_eo - matches[1].rm_so;
            result = (char*)malloc(length + 1);
            if (result) {
                strncpy(result, line + matches[1].rm_so, length);
                result[length] = '\0';
            }
        }
    }
    
    regfree(&regex);

    return result;
}

int _exists_in_subarray(int* arr, int start, int end, int value) {
    for (int i = start; i <= end; i++) {
        if (arr[i] == value)
            return 1;
        return 0;
    }
}

void delete_invidual(invidual* inv) {
    if (inv != NULL) {
        if (inv->gene.entries != NULL)
            free(inv->gene.entries);
        free(inv);
    }
}

problem* create_problem(char* file_path) {
    
    return "";
}

invidual* crossover_ox(invidual* inv1, invidual* inv2) {
    size_t inv_size = inv1->gene.size;
    
    int* child_inv = (int*)malloc(inv_size * sizeof(int));

    int start_sub_index = RANDOM(0, inv_size-1);
    int end_sub_index = RANDOM(start_sub_index, inv_size-1);

    size_t sub_size = end_sub_index - start_sub_index + 1;

    for (int i = 0; i < inv_size; i++) {
        child_inv[i] = BLANK;
    }

    for (int i = start_sub_index; i <= end_sub_index; i++) {
        child_inv[i] = inv1->gene.entries[i];
    }

    int current_index = 0;
    for (int i = 0; i < inv_size; i++) {
        int value = inv2->gene.entries[i];
        if (!_exists_in_subarray(child_inv, start_sub_index, end_sub_index, value)) {
            while (child_inv[current_index] != BLANK)
                current_index++;
            if (current_index < child_inv)
                child_inv[current_index] = value;
            else {
                return NULL; // error
            }
        }
    }

    return child_inv;
}

void mutation_swap(invidual* inv, double prob) {
    for (int i = 0; i < inv->gene.size; i++) {
        for (int j = 0; j < inv->gene.size; j++) {
            if (i == j) continue;
            if (rand() < prob) {
                int temp = inv->gene.entries[i];
                inv->gene.entries[i] = inv->gene.entries[j];
                inv->gene.entries[j] = temp;
            }
        }
    }
}