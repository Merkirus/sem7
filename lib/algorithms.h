#pragma once

#include "structure.h"
#include <limits.h>
#include "utils.h"
#include <float.h>
#include <stdbool.h>
#include <math.h>

invidual* selection_tournament(invidual** generation, size_t gen_size, size_t tournament_size);
void selection_elitism(invidual** generation, invidual** new_generation, size_t gen_size, size_t elite_size);
invidual* selection_roulette(invidual** generation, size_t gen_size, size_t placeholder);

invidual* crossover_ox(invidual* inv1, invidual* inv2);
invidual* crossover_pmx(invidual* inv1, invidual* inv2);
invidual* crossover_cx(invidual* inv1, invidual* inv2);

void mutation_swap(invidual* inv, double prob);
void mutation_single_swap(invidual* inv, int idx1, int idx2);
void mutation_inversion(invidual* inv, int idx1, int idx2);

void ant_traverse(ant* a, double** distance_matrix, double** total_pheromones, double alpha, double beta, size_t prob_size);
void update_pheromones(double** total_pheromones, ant* a, size_t prob_size);
void evaporate_pheromones(double** total_pheromones, double rho, size_t prob_size);

invidual** alg_greedy(problem* prob);
invidual** alg_random(problem* prob, int iter_size);
invidual*** alg_ea(problem* prob, int gen_size, int pop_size, int mutation_flag, int selection_flag, int elit_flag, invidual* (*selection)(invidual**,size_t,size_t), invidual* (*crossover)(invidual*,invidual*), void (*mutation)(invidual*,int,int), double px, double pm); 
invidual*** alg_ts(problem* prob, int iter_size, int neighbour_size, int tabu_size, int greedy_flag, void (*mutation)(invidual*,int,int));
invidual*** alg_sa(problem* prob, int iter_size, int* iter_count, int neighbour_size, int tabu_size, double start_temp, double end_temp, double alpha, void (*mutation)(invidual*,int,int));
invidual*** alg_aoc(problem* prob, int iter_size, int ant_size, int evaportation_flag, double alpha, double beta, double rho);
invidual*** alg_hybrid(problem* prob, int tabu_iter_size, int sa_iter_size, int tabu_neighbour_size, int sa_neighbour_size, int tabu_size, double start_temp, double end_temp, double alpha, int mutation_flag);