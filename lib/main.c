#include "algorithms.h"
#include "structure.h"
#include "utils.h"
#include <time.h>

#define GEN_SIZE 500
#define POP_SIZE 5000
#define ITERATIONS 35000
#define NEIGHBOUR_SIZE 5
#define TABU_SIZE 100

void run_ea(const char* problem_path, int gen_size, int pop_size, int elitism_flag, int selection_flag, int crossover_flag, int mutation_flag, double px, double pm) {
    srand(time(NULL));
    problem* prob = create_problem(problem_path);

    for (int i = 0; i < 10; i++) {
        invidual*** results;
        void (*mutation)(invidual*,int,int);
        invidual* (*crossover)(invidual*,invidual*);

        if (mutation_flag == 0)
            mutation = mutation_inversion;
        else
            mutation = mutation_single_swap;
        
        if (crossover_flag == 1)
            crossover = crossover_ox;
        else
            crossover = crossover_cx;
        
        results = alg_ea(prob, gen_size, pop_size, mutation_flag, selection_flag, elitism_flag, selection_tournament, crossover, mutation, px, pm);

        char buffer[16];
        sprintf(buffer, "output%d.csv", i);
        save_result(results, gen_size, pop_size, buffer);
    }
    
    delete_problem(prob);
}

void run_ts(const char* problem_path, int iter_size, int neighbour_size, int tabu_size, int mutation_flag, int greedy_flag) {
    srand(time(NULL));
    problem* prob = create_problem(problem_path);

    for (int i = 0; i < 10; i++) {
        invidual*** results;
        void (*mutation)(invidual*,int,int);

        if (mutation_flag == 0)
            mutation = mutation_inversion;
        else
            mutation = mutation_single_swap;
        
        results = alg_ts(prob, iter_size, neighbour_size, tabu_size, greedy_flag, mutation);
        char buffer[16];
        sprintf(buffer, "output%d.csv", i);
        save_result(results, iter_size, neighbour_size, buffer);
    }
    
    delete_problem(prob);

}

void run_sa(const char* problem_path, int iter_size, int neighbour_size, int tabu_size, double temp_start, double temp_end, double alpha, int mutation_flag) {
    srand(time(NULL));
    problem* prob = create_problem(problem_path);

    for (int i = 0; i < 10; i++) {
        int iter_count = 0;
        invidual*** results;
        void (*mutation)(invidual*,int,int);

        if (mutation_flag == 0)
            mutation = mutation_inversion;
        else
            mutation = mutation_single_swap;

        results = alg_sa(prob, iter_size, &iter_count, neighbour_size, tabu_size, temp_start, temp_end, alpha, mutation);

        char buffer[16];
        sprintf(buffer, "output%d.csv", i);
        save_result(results, 1, iter_count, buffer);
    }
    
    delete_problem(prob);
}

void run_aco(const char* problem_path, int iter_size, int ant_size, int evaporation_flag, double alpha, double beta, double rho) {
    srand(time(NULL));
    problem* prob = create_problem(problem_path);

    for (int i = 0; i < 10; i++) {
        invidual*** results;

        results = alg_aoc(prob, iter_size, ant_size, evaporation_flag, alpha, beta, rho);

        char buffer[16];
        sprintf(buffer, "output%d.csv", i);
        save_result(results, iter_size, ant_size, buffer);
    }

    printf("%d\n", EVAL_COUNT);
    
    delete_problem(prob);
}

void run_hybrid(const char* problem_path, int tabu_iter_size, int sa_iter_size, int tabu_neighbour_size, int sa_neighbour_size, int tabu_size, double start_temp, double end_temp, double alpha, int mutation_flag) {
    srand(time(NULL));
    problem* prob = create_problem(problem_path);


    for (int i = 0; i < 10; i++) {
        invidual*** results;
        void (*mutation)(invidual*,int,int);

        results = alg_hybrid(prob, tabu_iter_size, sa_iter_size, tabu_neighbour_size, sa_neighbour_size, tabu_size, start_temp, end_temp, alpha, mutation_flag);

        char buffer[16];
        sprintf(buffer, "output%d.csv", i);
        save_result(results, 1, tabu_iter_size, buffer);
    }

    printf("%d\n", EVAL_COUNT);
    
    delete_problem(prob);
}

int main() {
    // run_ea("../Prob/A-n48-k7.vrp", 1000, 50, 0, 5, 1, 0, 0.7, 0.2);
    // run_ts("../Prob/A-n48-k7.vrp", 10000, 20, 100, 0, 0);
    // run_sa("../Prob/A-n48-k7.vrp", 50000, 3, 100, 100.0, 0.001, 0.99995, 0);
    run_aco("../Prob/A-n48-k7.vrp", 500, 10, 1, 1.0, 5.0, 0.2);
    return 0;
}