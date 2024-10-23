#pragma once

#include <stdio.h>
#include <stdlib.h> 
#include <stddef.h>
#include <math.h>
#include <regex.h>
#include <string.h>
#include <limits.h>

#define RANDOM(min, max) (rand() % ((max) - (min) - (min) + 1) + (min))
#define BLANK -100
#define BUFFER 256
#define METADATA_SIZE 7
#define DEPOT 0

double _distance(double x1, double y1, double x2, double y2);
char* _get_file_data(char* line);
int _exists_in_subarray(int* arr, int start, int end, int value);
int _min_value_index(double* distances, size_t size);

typedef struct {
    size_t size;
    int* entries;
} vec2d;

typedef struct {
    size_t size;
    double* entries;
} d_vec2d;

typedef struct {
    size_t size_x;
    size_t size_y;
    int** entries;
} vec3d;

typedef struct {
    size_t size_x;
    size_t size_y;
    double** entries;
} d_vec3d;

typedef struct {
    vec2d gene;
} invidual;

typedef struct {
    int cap;
    d_vec3d distance_matrix;
    vec2d demand;
} problem;


invidual* create_invidual_empty(size_t gene_size);
invidual* create_invidual_random(size_t gene_size);
invidual* create_invidual_greedy(int start, int* closest_stations, size_t gene_size);
void delete_invidual(invidual* inv);
double evaluate_invidual(problem* prob, invidual* inv);

problem* create_problem(char* file_path);

invidual* crossover_ox(invidual* inv1, invidual* inv2);
invidual* crossover_pmx(invidual* inv1, invidual* inv2);
invidual* crossover_cx(invidual* inv1, invidual* inv2);

void mutation_swap(invidual* inv, double prob);

double alg_greedy(problem* prob);
double alg_random(problem* prob);
double alg_ea(problem* prob);

double _distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
}

char* _get_file_data(char* line) {
    regex_t regex;
    regmatch_t matches[2];
    char* pattern = "\\w+\\s*:\\s*(.*)";
    char* result = NULL;
    
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

int _min_value_index(double* distances, size_t size) {
    int min = INT_MAX;
    int min_index = 0;

    for (int i = 0; i < size; i++) {
        int value = distances[i];

        if (value != 0 && value < min) {
            min = value;
            min_index = i;
        }
    }

    return min_index;
}

invidual* create_invidual_empty(size_t gene_size) {
    return NULL;
}

invidual* create_invidual_random(size_t gene_size) {
    return NULL;
}

invidual* create_invidual_greedy(int start, int* closest_stations, size_t gene_size) {
    return NULL;
}

void delete_invidual(invidual* inv) {
    if (inv != NULL) {
        if (inv->gene.entries != NULL)
            free(inv->gene.entries);
        free(inv);
    }
}

double evaluate_invidual(problem* prob, invidual* inv) {
    int curr_cap = prob->cap;
    int curr_station = inv->gene.entries[0];

    double distance = 0.0;
    int route[BUFFER];
    for (int i = 0; i < BUFFER; i++) route[i] = BLANK;

    for (int i = 1, route_index = 0; i < inv->gene.size; i++) {
        int next_station = inv->gene.entries[i];
        int curr_demand = prob->demand.entries[next_station];
        if (curr_cap - curr_demand < 0) {
            distance += prob->distance_matrix.entries[curr_station][DEPOT];
            curr_cap = prob->cap;
            curr_station = DEPOT;
            i--;
            route[route_index] = DEPOT;
            route_index++;
            continue;
        }
        curr_cap -= curr_demand;
        distance += prob->distance_matrix.entries[curr_station][next_station];
        route[route_index] = next_station;
        route_index++;
        curr_station = next_station;
    }

    distance += prob->distance_matrix.entries[curr_station][DEPOT];

    return distance;
}

problem* create_problem(char* file_path) {
    problem* prob = (problem*)malloc(sizeof(problem));
    FILE* file;
    char line[BUFFER];

    file = fopen(file_path, "r");

    if (file == NULL)
        return NULL;

    // filetering junk data
    for (int i = 0; i < METADATA_SIZE; i++) {
        if (!fgets(line, sizeof(line), file))
            return NULL;

        char* value = _get_file_data(line);
        
        if (value != NULL) {
            int i_value;
            if (i == 3) {
                i_value = atoi(value);
                prob->distance_matrix.size_x = i_value;
                prob->distance_matrix.size_y = i_value;
                prob->demand.size = i_value;
            }
            else if (i == 5)
                i_value = atoi(value);
                prob->cap = i_value;
        }
    }

    prob->distance_matrix.entries[prob->distance_matrix.size_x][prob->distance_matrix.size_y];
    for (int i = 0; i < prob->distance_matrix.size_x; i++)
        prob->distance_matrix.entries[i] = (double*)malloc(prob->distance_matrix.size_y * sizeof(double));

    prob->demand.entries = (int*)malloc(prob->demand.size * sizeof(int));
    
    int coordinates[prob->distance_matrix.size_x][2];

    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        if (!fgets(line, sizeof(line), file))
            return NULL;

        int* _;
        if (!sscanf(line, "%d %d %d", _, &coordinates[i][0], &coordinates[i][1]))
            return NULL;
    }

    if (!fgets(line, sizeof(line), file)) // Demand section
        return NULL;

    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        if (!fgets(line, sizeof(line), file))
            return NULL;
        
        int* _;
        if (!sscanf(line, "%d %d", _, &prob->demand.entries[i]))
            return NULL;
    }

    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        for (int j = 0; j < prob->distance_matrix.size_y; j++) {
            int x1 = coordinates[i][0];
            int y1 = coordinates[i][1];
            int x2 = coordinates[j][0];
            int y2 = coordinates[j][1];
            prob->distance_matrix.entries[i][j] = _distance(x1,y1,x2,y2);
        }
    }

    fclose(file);

    return prob;
}

invidual* crossover_ox(invidual* inv1, invidual* inv2) {
    size_t inv_size = inv1->gene.size;
    
    invidual* child_inv = create_invidual_empty(inv_size);

    int start_sub_index = RANDOM(0, inv_size-1);
    int end_sub_index = RANDOM(start_sub_index, inv_size-1);

    size_t sub_size = end_sub_index - start_sub_index + 1;

    for (int i = 0; i < inv_size; i++) {
        child_inv->gene.entries[i] = BLANK;
    }

    for (int i = start_sub_index; i <= end_sub_index; i++) {
        child_inv->gene.entries[i] = inv1->gene.entries[i];
    }

    int current_index = 0;
    for (int i = 0; i < inv_size; i++) {
        int value = inv2->gene.entries[i];
        if (!_exists_in_subarray(child_inv->gene.entries, start_sub_index, end_sub_index, value)) {
            while (child_inv->gene.entries[current_index] != BLANK)
                current_index++;
            if (current_index < child_inv->gene.entries)
                child_inv->gene.entries[current_index] = value;
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

double alg_greedy(problem* prob) {
    int closest_stations[prob->distance_matrix.size_x];

    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        closest_stations[i] = _min_value_index(prob->distance_matrix.entries[i], prob->distance_matrix.size_y);
    }

    invidual** invs = (invidual**)malloc(prob->distance_matrix.size_x * sizeof(invidual*));
    for (int i = 0; i < prob->distance_matrix.size_x; i++) invs[i] = create_invidual_greedy(i, closest_stations, prob->distance_matrix.size_x);

    double costs[prob->distance_matrix.size_x];
    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        costs[i] = evaluate_invidual(prob, invs[i]);
    }

    for (int i = 0; i < prob->distance_matrix.size_x; i++)
        delete_invidual(invs[i]);

    int result = costs[_min_value_index(costs, prob->distance_matrix.size_x)];
    
    free(invs);

    return result;
}

double alg_random(problem* prob) {
    invidual* inv = create_invidual_random(prob->distance_matrix.size_x);
    int result = evaluate_invidual(prob, inv);
    delete_invidual(inv);
    return result;
}

double alg_ea(problem* prob) {

    return 0.0;
}