#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <regex.h>

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
    double fitness;
} invidual;

typedef struct {
    int cap;
    d_vec3d distance_matrix;
    vec2d demand;
} problem;

typedef struct {
    invidual* solution;
    d_vec3d pheromones;
    vec2d visited;
    int current_position;
    int current_index;
    int starting_position;
} ant;

typedef struct {
    int move_type;
    int idx1;
    int idx2;
    int lifespan;
} tabu_move;

typedef struct {
    tabu_move** moves;
    size_t moves_size;
    int current_size;
} tabu_list;

extern int EVAL_COUNT;

invidual* create_invidual_empty(size_t gene_size);
invidual* create_invidual_random(size_t gene_size);
invidual* create_invidual_greedy(int start, int** closest_stations, size_t gene_size);
invidual* copy_invidual(invidual* inv);
void delete_invidual(invidual* inv);
double evaluate_invidual(problem* prob, invidual* inv);

ant* create_ant(int starting_position, size_t gene_size);
void delete_ant(ant* a);
void move_ant(ant* a, int destination);

tabu_move* create_tabu_move(int move_type, int idx1, int idx2, int lifespan);
tabu_move* copy_tabu_move(tabu_move* tm);
tabu_list* create_tabu_list(size_t tabu_size);
void add_tabu_move(tabu_list* tl, tabu_move* tm);
int is_tabu(tabu_list* tl, tabu_move* tm);
void update_tabu_list(tabu_list* tl);
void delete_tabu_list(tabu_list* tl);

problem* create_problem(char* file_path);
void delete_problem(problem* problem);