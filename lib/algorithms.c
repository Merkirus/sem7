#include "algorithms.h"
#include "time.h"

#define RANDOM(min, max) (min + rand() % ((max) - (min) + 1))
#define PROB() ((double)rand() / RAND_MAX)
#define BLANK -100
#define BUFFER 256
#define BUFFER_SA 500000
#define METADATA_SIZE 7
#define DEPOT 0
#define TOUR_SIZE 100
#define ELITE_SIZE 50
#define PX 0.3
#define PM 0.15
#define TABU_TENURE 10

invidual* selection_tournament(invidual** generation, size_t gen_size, size_t tournament_size) {
    int best_index = RANDOM(0, gen_size - 1);
    double best_fitness = generation[best_index]->fitness;

    for (int i = 1; i < tournament_size; i++) {
        int random_index = RANDOM(0, gen_size - 1);
        double current_fitness = generation[random_index]->fitness;

        if (current_fitness < best_fitness) {
            best_fitness = current_fitness;
            best_index = random_index;
        }
    }

    return generation[best_index];
}

void selection_elitism(invidual** generation, invidual** new_generation, size_t gen_size, size_t elite_size) {
    qsort(generation, gen_size, sizeof(invidual*), compare_fitness);

    for (int i = 0; i < elite_size; i++) {
        new_generation[i] = generation[i];
    }
}

invidual* selection_roulette(invidual** generation, size_t gen_size, size_t placeholder) {
    // Calculate the total fitness of the population
    double total_fitness = 0.0;
    for (int i = 0; i < gen_size; i++) {
        total_fitness += generation[i]->fitness;
    }

    // Generate a random number between 0 and total_fitness
    double rand_value = ((double)rand() / RAND_MAX) * total_fitness;

    // Select an individual based on the random value
    double accumulated_fitness = 0.0;
    for (int i = 0; i < gen_size; i++) {
        accumulated_fitness += generation[i]->fitness;
        if (accumulated_fitness >= rand_value) {
            return generation[i];
        }
    }

    return generation[gen_size - 1];  // Fallback (should never reach here)
}

// invidual* crossover_ox(invidual* inv1, invidual* inv2) {
//     size_t inv_size = inv1->gene.size;
    
//     invidual* child_inv = create_invidual_empty(inv_size);

//     int start_sub_index = RANDOM(0, inv_size - 1);
//     int end_sub_index = RANDOM(start_sub_index, inv_size - 1);

//     for (int i = 0; i < inv_size; i++) {
//         child_inv->gene.entries[i] = BLANK;
//     }

//     for (int i = start_sub_index; i <= end_sub_index; i++) {
//         child_inv->gene.entries[i] = inv1->gene.entries[i];
//     }

//     int current_index = 0;
//     for (int i = 0; i < inv_size; i++) {
//         int value = inv2->gene.entries[i];
//         int gene_present = 0;
//         for (int j = start_sub_index; j <= end_sub_index; j++) {
//             if (inv2->gene.entries[i] == child_inv->gene.entries[j]) {
//                 gene_present = 1;
//                 break;
//             }
//         }

//         if (!gene_present) {
//             while (current_index >= start_sub_index && current_index <= end_sub_index)
//                 current_index++;

//             child_inv->gene.entries[current_index++] = inv2->gene.entries[i];
//         }
//     }

//     return child_inv;
// }

invidual* crossover_ox(invidual* inv1, invidual* inv2) {
    size_t inv_size = inv1->gene.size;
    
    invidual* child_inv = create_invidual_empty(inv_size);

    int start_sub_index = RANDOM(0, inv_size - 1);
    int end_sub_index = RANDOM(start_sub_index, inv_size - 1);

    for (int i = 0; i < inv_size; i++) child_inv->gene.entries[i] = BLANK;
    int* gene_present = (int*)calloc(inv_size, sizeof(int));

    for (int i = start_sub_index; i <= end_sub_index; i++) {
        int value = inv1->gene.entries[i];
        child_inv->gene.entries[i] = value;
        gene_present[value] = 1;
    }

    int current_index = 0;
    for (int i = 0; i < inv_size; i++) {
        int value = inv2->gene.entries[i];

        if (!gene_present[value]) {
            while (child_inv->gene.entries[current_index] != BLANK) {
                current_index++;
            }
            child_inv->gene.entries[current_index++] = value;
        }
    }

    free(gene_present);

    return child_inv;
}

invidual* crossover_pmx(invidual* inv1, invidual* inv2) {
    size_t inv_size = inv1->gene.size;
    
    invidual* child_inv = create_invidual_empty(inv_size);

    int start_sub_index = RANDOM(0, inv_size-1);
    int end_sub_index = RANDOM(start_sub_index, inv_size-1);

    for (int i = 0; i < inv_size; i++)
        child_inv->gene.entries[i] = BLANK;

    for (int i = start_sub_index; i <= end_sub_index; i++)
        child_inv->gene.entries[i] = inv1->gene.entries[1];

    for (int i = 0; i < inv2->gene.size; i++) {
        int gene_present = 0;
        for (int j = start_sub_index; j <= end_sub_index; j++) {
            if (inv2->gene.entries[i] == child_inv->gene.entries[j]) {
                gene_present = 1;
                break;
            }
        }

        if (!gene_present) {
            int idx = -1;
            for (int j = 0; j < inv1->gene.size; j++) {
                if (child_inv->gene.entries[j] == BLANK) {
                    idx = j;
                    break;
                }
            }

            child_inv->gene.entries[idx] = inv2->gene.entries[i];
        }
    }

    return child_inv;
}

invidual* crossover_cx(invidual* inv1, invidual* inv2) {
    int size = inv1->gene.size;
    invidual *child = (invidual*)malloc(sizeof(invidual));
    child->gene.entries = (int*)malloc(size * sizeof(int));
    child->gene.size = size;
    
    // Initialize child genes with a special marker (-1) to indicate unfilled positions
    for (int i = 0; i < size; i++) {
        child->gene.entries[i] = -1;
    }

    bool *visited = (bool*)calloc(size, sizeof(bool)); // Track visited indices in cycles
    int cycle_start = 0; // Start with the first gene position
    bool use_inv1 = true; // Alternate cycles between parents

    while (cycle_start < size) {
        if (visited[cycle_start]) {
            cycle_start++;
            continue;
        }

        int pos = cycle_start;
        
        // Perform the cycle, marking visited genes and copying from one parent
        do {
            visited[pos] = true;
            child->gene.entries[pos] = use_inv1 ? inv1->gene.entries[pos] : inv2->gene.entries[pos];
            
            // Find the index in inv1 where inv2[pos] is located (to continue the cycle)
            int gene = inv2->gene.entries[pos];
            for (int i = 0; i < size; i++) {
                if (inv1->gene.entries[i] == gene) {
                    pos = i;
                    break;
                }
            }
        } while (pos != cycle_start); // Complete the cycle when back at the start

        // Alternate to the other parent for the next cycle
        use_inv1 = !use_inv1;
        cycle_start++;
    }

    free(visited);
    return child;
}

void mutation_swap(invidual* inv, double prob) {
    for (int i = 0; i < inv->gene.size; i++) {
        if (PROB() < prob) {
            int j = RANDOM(0, inv->gene.size - 1);
            int k = RANDOM(0, inv->gene.size - 1);

            while (j == k)
                k = RANDOM(0, inv->gene.size - 1);
            
            int temp = inv->gene.entries[j];
            inv->gene.entries[j] = inv->gene.entries[k];
            inv->gene.entries[k] = temp;
        }
    }
}

void mutation_single_swap(invidual* inv, int idx1, int idx2) {
    int temp = inv->gene.entries[idx1];
    inv->gene.entries[idx1] = inv->gene.entries[idx2];
    inv->gene.entries[idx2] = temp;
}

void mutation_inversion(invidual* inv, int idx1, int idx2) {
    size_t inv_size = inv->gene.size;
    int start_sub_index = RANDOM(0, inv_size-1);
    int end_sub_index = RANDOM(start_sub_index, inv_size-1);
    _reverse_subarray(inv->gene.entries, idx1, idx2);
}

void update_pheromones(double** total_pheromones, ant* a, size_t prob_size) {
    double pheromone_value = 10000.0 / a->solution->fitness;
    for (int i = 0; i < prob_size - 1; i++) {
        int from = a->solution->gene.entries[i];
        int to = a->solution->gene.entries[i+1];
        total_pheromones[from][to] += pheromone_value;
        total_pheromones[to][from] += pheromone_value;
    }

    int from = a->solution->gene.entries[prob_size - 1];
    int to = a->solution->gene.entries[0];
    total_pheromones[from][to] += pheromone_value;
    total_pheromones[to][from] += pheromone_value;
}

void evaporate_pheromones(double** total_pheromones, double rho, size_t prob_size) {
    for (int i = 0; i < prob_size; i++)
        for (int j = 0; j < prob_size; j++)
        total_pheromones[i][j] *= (1 - rho);
    
}

void ant_traverse(ant* a, double** distance_matrix, double** total_pheromones, double alpha, double beta, size_t prob_size) {
    for (int i = 1; i < prob_size; i++) { // Loop through all steps except the starting point
        double traversal[prob_size];
        double total_traversal = 0.0;

        // Calculate probabilities for all possible moves
        for (int j = 0; j < prob_size; j++) {
            if (a->visited.entries[j]) {
                traversal[j] = 0; // Already visited nodes are not eligible
            } else {
                double distance = distance_matrix[a->current_position][j];
                double pheromone = total_pheromones[a->current_position][j];

                // Ensure valid inputs
                if (distance <= 0.0) distance = 1e-10; // Prevent divide by zero
                if (pheromone <= 0.0) pheromone = 1e-10; // Prevent invalid pow()

                double v1 = pow(pheromone, alpha);
                double v2 = pow(1.0 / distance, beta);
                traversal[j] = v1 * v2;
                total_traversal += traversal[j];
            }
        }

        // Normalize probabilities
        if (total_traversal > 0) {
            for (int j = 0; j < prob_size; j++) {
                traversal[j] /= total_traversal;
            }
        } else {
            fprintf(stderr, "Error: Total traversal probability is zero!\n");
            exit(EXIT_FAILURE);
        }

        // Select destination based on probability
        double random_value = PROB(); // Generate a random number between 0 and 1
        double cumulative_probability = 0.0;
        int dest = -1;

        for (int j = 0; j < prob_size; j++) {
            if (!a->visited.entries[j]) {
                cumulative_probability += traversal[j];
                if (random_value <= cumulative_probability) {
                    dest = j;
                    break;
                }
            }
        }

        // Fallback: If no destination was chosen (should not happen)
        if (dest == -1) {
            fprintf(stderr, "Warning: No valid destination selected using probabilities. Using fallback.\n");
            double max_probability = -1.0;
            for (int j = 0; j < prob_size; j++) {
                if (!a->visited.entries[j] && traversal[j] > max_probability) {
                    max_probability = traversal[j];
                    dest = j;
                }
            }
        }

        // Move the ant to the chosen destination
        if (dest != -1) {
            move_ant(a, dest);
        } else {
            fprintf(stderr, "Error: No valid move available for the ant.\n");
            exit(EXIT_FAILURE);
        }
    }
}

// void ant_traverse(ant* a, double** distance_matrix, double** total_pheromones, double alpha, double beta, size_t prob_size) {
//     for (int i = 1; i < prob_size; i++) { // - 1 because we skip starting_point
//         double traversal[prob_size];
//         double total_travesal = 0.0;

//         for (int j = 0; j < prob_size; j++) {
//             if (a->visited.entries[j])
//                 traversal[j] = 0;
//             else {
//                 double distance = distance_matrix[a->current_position][j];
//                 double pheromone = total_pheromones[a->current_position][j];
//                 double v1 = pow(pheromone, alpha);
//                 double v2 = pow(1.0 / distance, beta);
//                 double path_value = v1 * v2;
//                 traversal[j] = path_value;
//                 total_travesal += path_value;
//             }
//         }

//         for (int j = 0; j < prob_size; j++) {
//             traversal[j] /= total_travesal; // normalizing
//         }

//         int dest = -1;

//         double add_traversal = 0.0;
//         for (int j = 0; j < prob_size; j++) {
//             if (!a->visited.entries[j]) {
//                 add_traversal += traversal[j];
//                 if (add_traversal >= PROB()) {
//                     dest = j;
//                     break;
//                 }
//             }
//         }

//         if (dest == -1) {
//             double max_value = -1.0;
//             int index = 0;
//             for (int j = 0; j < prob_size; j++) {
//                 if (!a->visited.entries[j]) {
//                     if (max_value < traversal[j]) {
//                         max_value = traversal[j];
//                         index = j;
//                     }
//                 }
//             }
//             dest = index;
//         }

//         move_ant(a, dest);
//     }
// }

invidual** alg_greedy(problem* prob) {
    size_t size = prob->distance_matrix.size_x;
    
    invidual** results = (invidual**)malloc(size * sizeof(invidual*));
    for (int i = 0; i < size; i++)
        results[i] = (invidual*)malloc(sizeof(invidual));
    
    double temp_arr[size][size];
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            temp_arr[i][j] = prob->distance_matrix.entries[i][j];
        
    for (int i = 0; i < size; i++)
        qsort(temp_arr[i], size, sizeof(double), compare);

    int** closest_stations = (int**)malloc(size * sizeof(int*));
    for (int i = 0; i < size; i++)
        closest_stations[i] = (int*)malloc(size * sizeof(int));

    for (int m = 0; m < size; m++)
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                if (temp_arr[m][i] == prob->distance_matrix.entries[m][j])
                    closest_stations[m][i] = j;

    for (int i = 0; i < size; i++) {
        invidual* inv = create_invidual_greedy(i, closest_stations, size);
        inv->fitness = evaluate_invidual(prob, inv);
        results[i] = inv;
    }

    for (int i = 0; i < size; i++)
        free(closest_stations[i]);

    free(closest_stations);

    return results;
}

invidual** alg_random(problem* prob, int iter_size) {
    invidual** results = (double*)malloc(iter_size * sizeof(double));

    for (int i = 0; i < iter_size; i++) {
        invidual* inv = create_invidual_random(prob->distance_matrix.size_x);
        inv->fitness = evaluate_invidual(prob, inv);
        results[i] = inv;
    }

    return results;
}

invidual*** alg_ea(problem* prob, int gen_size, int pop_size, int mutation_flag, int selection_flag, int elit_flag, invidual* (*selection)(invidual**,size_t,size_t), invidual* (*crossover)(invidual*,invidual*), void (*mutation)(invidual*,int,int), double px, double pm) {
    invidual*** generations = (invidual***)malloc(gen_size * sizeof(invidual**));
    for (int i = 0; i < gen_size; i++)
        generations[i] = (invidual**)malloc(pop_size * sizeof(invidual*));

    for (int i = 0; i < gen_size; i++)
        for (int j = 0; j < pop_size; j++)
            generations[i][j] = (invidual*)malloc(sizeof(invidual));

    for (int i = 0; i < pop_size; i++)
        generations[0][i] = create_invidual_random(prob->distance_matrix.size_x);

    for (int i = 0; i < pop_size; i++) {
        int result = evaluate_invidual(prob, generations[0][i]);
        generations[0][i]->fitness = result;
    }

    for (int gen_iter = 0; gen_iter < gen_size-1; gen_iter++) {
        
        if (elit_flag != 0)
            selection_elitism(generations[gen_iter], generations[gen_iter+1], pop_size, elit_flag);

        for (int pop_iter = elit_flag; pop_iter < pop_size; pop_iter++) {

            invidual* inv1;
            invidual* inv2;

            inv1 = selection(generations[gen_iter], pop_size, selection_flag);
            inv2 = selection(generations[gen_iter], pop_size, selection_flag);

            invidual* child_inv;

            if (PROB() < px)
                child_inv = crossover(inv1, inv2);
            else
                child_inv = inv1;

            if (mutation_flag == 1)
                mutation_swap(child_inv, pm);
            else {
                if (PROB() < pm) {
                    int idx1 = RANDOM(0, child_inv->gene.size - 1);
                    int idx2 = RANDOM(idx1, child_inv->gene.size - 1);
                    mutation(child_inv, idx1, idx2);
                }
            }

            int result = evaluate_invidual(prob, child_inv);
            child_inv->fitness = result;

            generations[gen_iter+1][pop_iter] = child_inv;
        }
    }

    return generations;
}

// invidual*** alg_ea(problem* prob, int gen_size, int pop_size, int mutation_flag, int selection_flag, int elit_flag, invidual* (*selection)(invidual**,size_t,size_t), invidual* (*crossover)(invidual*,invidual*), void (*mutation)(invidual*,int,int), double px, double pm) {
//     invidual*** generations = (invidual***)malloc(gen_size * sizeof(invidual**));
//     for (int i = 0; i < gen_size; i++)
//         generations[i] = (invidual**)malloc(pop_size * sizeof(invidual*));

//     hash_table* hs = create_hash_table();

//     for (int i = 0; i < gen_size; i++)
//         for (int j = 0; j < pop_size; j++)
//             generations[i][j] = (invidual*)malloc(sizeof(invidual));

//     for (int i = 0; i < pop_size; i++)
//         generations[0][i] = create_invidual_random(prob->distance_matrix.size_x);

//     for (int i = 0; i < pop_size; i++) {
//         double is_fitness = search(hs, generations[0][i]->gene.entries, generations[0][i]->gene.size);
//         if (isnan(is_fitness)) {
//             int result = evaluate_invidual(prob, generations[0][i]);
//             generations[0][i]->fitness = result;
//             insert(hs, generations[0][i]->gene.entries, generations[0][i]->gene.size, result);
//         } else
//             generations[0][i]->fitness = is_fitness;
//     }

//     for (int gen_iter = 0; gen_iter < gen_size-1; gen_iter++) {
        
//         if (elit_flag != 0)
//             selection_elitism(generations[gen_iter], generations[gen_iter+1], pop_size, elit_flag);

//         for (int pop_iter = elit_flag; pop_iter < pop_size; pop_iter++) {

//             invidual* inv1;
//             invidual* inv2;

//             inv1 = selection(generations[gen_iter], pop_size, selection_flag);
//             inv2 = selection(generations[gen_iter], pop_size, selection_flag);

//             invidual* child_inv;

//             if (PROB() < px)
//                 child_inv = crossover(inv1, inv2);
//             else
//                 child_inv = inv1;

//             if (mutation_flag == 1)
//                 mutation_swap(child_inv, pm);
//             else {
//                 if (PROB() < pm) {
//                     int idx1 = RANDOM(0, child_inv->gene.size - 1);
//                     int idx2 = RANDOM(idx1, child_inv->gene.size - 1);
//                     mutation(child_inv, idx1, idx2);
//                 }
//             }

//             double is_fitness = search(hs, child_inv->gene.entries, child_inv->gene.size);
//             if (isnan(is_fitness)) {
//                 int result = evaluate_invidual(prob, child_inv);
//                 child_inv->fitness = result;
//                 insert(hs, child_inv->gene.entries, child_inv->gene.size, result);
//             } else
//                 child_inv->fitness = is_fitness;

//             generations[gen_iter+1][pop_iter] = child_inv;
//         }
//     }

//     return generations;
// }

invidual*** alg_ts(problem* prob, int iter_size, int neighbour_size, int tabu_size, int greedy_flag, void (*mutation)(invidual*,int,int)) {
    invidual* starting_inv;
    if (greedy_flag == 1) {
        size_t size = prob->distance_matrix.size_x;

        double temp_arr[size][size];
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                temp_arr[i][j] = prob->distance_matrix.entries[i][j];
            
        for (int i = 0; i < size; i++)
            qsort(temp_arr[i], size, sizeof(double), compare);

        int** closest_stations = (int**)malloc(size * sizeof(int*));
        for (int i = 0; i < size; i++)
            closest_stations[i] = (int*)malloc(size * sizeof(int));

        for (int m = 0; m < size; m++)
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    if (temp_arr[m][i] == prob->distance_matrix.entries[m][j])
                        closest_stations[m][i] = j;

        starting_inv = create_invidual_greedy(DEPOT, closest_stations, size);

        for (int i = 0; i < size; i++)
            free(closest_stations[i]);

        free(closest_stations);
    }
    else
        starting_inv = create_invidual_random(prob->distance_matrix.size_x);

    starting_inv->fitness = evaluate_invidual(prob, starting_inv);

    invidual* current_inv = copy_invidual(starting_inv);
    // invidual* best_solution;
    // best_solution = copy_invidual(current_solution);

    invidual*** results = (invidual***)malloc(iter_size * sizeof(invidual**));

    for (int i = 0; i < iter_size; i++)
        results[i] = (invidual**)malloc(neighbour_size * sizeof(invidual*));

    for (int i = 0; i < iter_size; i++)
        for (int j = 0; j < neighbour_size; j++)
            results[i][j] = (invidual*)malloc(sizeof(invidual));

    tabu_list* global_tabu_list = create_tabu_list(tabu_size);
    tabu_list* local_tabu_list = create_tabu_list(neighbour_size);

    for (int i = 0; i < iter_size; i++) {
        if (local_tabu_list == NULL)
            local_tabu_list = create_tabu_list(neighbour_size);

        invidual** neighbours = (invidual**)malloc(neighbour_size * sizeof(invidual*));
        for (int j = 0; j < neighbour_size; j++) neighbours[j] = (invidual*)malloc(sizeof(invidual));
        int neighbour_count = 0;
        int best_move_idx = 0;

        while (neighbour_count < neighbour_size) {
            int idx1 = RANDOM(0, current_inv->gene.size - 1);
            int idx2 = RANDOM(idx1, current_inv->gene.size - 1);
            tabu_move* tm = create_tabu_move(0, idx1, idx2, TABU_TENURE);
            if (!is_tabu(global_tabu_list, tm) || !is_tabu(local_tabu_list, tm)) {
                neighbours[neighbour_count] = copy_invidual(current_inv);
                mutation(neighbours[neighbour_count], idx1, idx2);
                neighbours[neighbour_count]->fitness = evaluate_invidual(prob, neighbours[neighbour_count]);
                neighbour_count++;
                add_tabu_move(local_tabu_list, tm);
            }
        }

        invidual* best_neighbour = NULL;

        for (int j = 0; j < neighbour_count; j++) {
            if (best_neighbour == NULL || neighbours[j]->fitness < best_neighbour->fitness) {
                best_neighbour = neighbours[j];
                best_move_idx = j;
            }
        }

        if (current_inv != NULL) {
            delete_invidual(current_inv);
            current_inv = NULL;
        }

        current_inv = copy_invidual(best_neighbour);
        // delete_invidual(current_solution);
        // current_solution = copy_invidual(best_neighbour);

        tabu_move* best_move = copy_tabu_move(local_tabu_list->moves[best_move_idx]);

        add_tabu_move(global_tabu_list, best_move);

        update_tabu_list(global_tabu_list);

        results[i] = neighbours;

        delete_tabu_list(local_tabu_list);
        local_tabu_list = NULL;
    }

    delete_invidual(starting_inv);

    delete_tabu_list(global_tabu_list);
    global_tabu_list = NULL;

    return results;
}

invidual*** alg_sa(problem* prob, int iter_size, int* iter_count, int neighbour_size, int tabu_size, double start_temp, double end_temp, double alpha, void (*mutation)(invidual*,int,int)) {
    invidual* current_inv = create_invidual_random(prob->distance_matrix.size_x);
    current_inv->fitness = evaluate_invidual(prob, current_inv);

    double temp = start_temp;

    invidual*** results = (invidual***)malloc(sizeof(invidual**));
    results[0] = (invidual**)malloc(iter_size * sizeof(invidual*));
    for (int i = 0; i < iter_size; i++)
        results[0][i] = (invidual*)malloc(sizeof(invidual));
    
    tabu_list* global_tabu_list = create_tabu_list(tabu_size);
    tabu_list* local_tabu_list = create_tabu_list(tabu_size);
    
    while (temp > end_temp && *iter_count < iter_size) {
        if (local_tabu_list == NULL)
            local_tabu_list = create_tabu_list(tabu_size);

        invidual** neighbours = (invidual**)malloc(neighbour_size * sizeof(invidual*));
        for (int i = 0; i < neighbour_size; i++) neighbours[i] = (invidual*)malloc(sizeof(invidual));
        int neighbour_count = 0;
        int best_move_idx = 0;

        while (neighbour_count < neighbour_size) {
            int idx1 = RANDOM(0, current_inv->gene.size - 1);
            int idx2 = RANDOM(idx1, current_inv->gene.size - 1);
            tabu_move* tm = create_tabu_move(0, idx1, idx2, TABU_TENURE);
            if (!is_tabu(global_tabu_list, tm) || !is_tabu(local_tabu_list, tm)) {
                neighbours[neighbour_count] = copy_invidual(current_inv);
                mutation(neighbours[neighbour_count], idx1, idx2);
                neighbours[neighbour_count]->fitness = evaluate_invidual(prob, neighbours[neighbour_count]);
                neighbour_count++;
                add_tabu_move(local_tabu_list, tm);
            }
        }

        invidual* best_neighbour = NULL;

        for (int j = 0; j < neighbour_count; j++) {
            if (best_neighbour == NULL || neighbours[j]->fitness < best_neighbour->fitness) {
                best_neighbour = neighbours[j];
                best_move_idx = j;
            }
        }

        double delta = best_neighbour->fitness - current_inv->fitness;
        
        if (delta < 0 || exp(-delta / temp) > PROB()) {
            current_inv = NULL;

            current_inv = copy_invidual(best_neighbour);

            results[0][*iter_count] = current_inv;

            (*iter_count)++;

            tabu_move* best_move = copy_tabu_move(local_tabu_list->moves[best_move_idx]);

            add_tabu_move(global_tabu_list, best_move);

            update_tabu_list(global_tabu_list);
        }

        delete_tabu_list(local_tabu_list);
        local_tabu_list = NULL;

        for (int i = 0; i < neighbour_size; i++)
            delete_invidual(neighbours[i]);

        free(neighbours);

        temp *= alpha;
    } 

    return results;
}

// invidual*** alg_sa(problem* prob, int iter_size, int* iter_count, int neighbour_size, int tabu_size, double start_temp, double end_temp, double alpha, void (*mutation)(invidual*,int,int)) {
//     invidual* current_inv = create_invidual_random(prob->distance_matrix.size_x);
//     hash_table* hs = create_hash_table();
//     current_inv->fitness = evaluate_invidual(prob, current_inv);

//     insert(hs, current_inv->gene.entries, current_inv->gene.size, current_inv->fitness);

//     double temp = start_temp;

    
//     invidual*** results = (invidual***)malloc(sizeof(invidual**));
//     results[0] = (invidual**)malloc(iter_size * sizeof(invidual*));
//     for (int i = 0; i < iter_size; i++)
//         results[0][i] = (invidual*)malloc(sizeof(invidual));
    
//     tabu_list* global_tabu_list = create_tabu_list(tabu_size);
//     tabu_list* local_tabu_list = create_tabu_list(tabu_size);
    
//     while (temp > end_temp && *iter_count < iter_size) {
//         if (local_tabu_list == NULL)
//             local_tabu_list = create_tabu_list(tabu_size);

//         invidual** neighbours = (invidual**)malloc(neighbour_size * sizeof(invidual*));
//         for (int i = 0; i < neighbour_size; i++) neighbours[i] = (invidual*)malloc(sizeof(invidual));
//         int neighbour_count = 0;
//         int best_move_idx = 0;

//         while (neighbour_count < neighbour_size) {
//             int idx1 = RANDOM(0, current_inv->gene.size - 1);
//             int idx2 = RANDOM(idx1, current_inv->gene.size - 1);
//             tabu_move* tm = create_tabu_move(0, idx1, idx2, TABU_TENURE);
//             if (!is_tabu(global_tabu_list, tm) || !is_tabu(local_tabu_list, tm)) {
//                 neighbours[neighbour_count] = copy_invidual(current_inv);
//                 mutation(neighbours[neighbour_count], idx1, idx2);
//                 double is_fitness = search(hs, neighbours[neighbour_count]->gene.entries, neighbours[neighbour_count]->gene.size);
//                 if (isnan(is_fitness)) {
//                     neighbours[neighbour_count]->fitness = evaluate_invidual(prob, neighbours[neighbour_count]);
//                     insert(hs, neighbours[neighbour_count]->gene.entries, neighbours[neighbour_count]->gene.size, neighbours[neighbour_count]->fitness);
//                 } else
//                     neighbours[neighbour_count]->fitness = is_fitness;
//                 neighbour_count++;
//                 add_tabu_move(local_tabu_list, tm);
//             }
//         }

//         invidual* best_neighbour = NULL;

//         for (int j = 0; j < neighbour_count; j++) {
//             if (best_neighbour == NULL || neighbours[j]->fitness < best_neighbour->fitness) {
//                 best_neighbour = neighbours[j];
//                 best_move_idx = j;
//             }
//         }

//         double delta = best_neighbour->fitness - current_inv->fitness;
        
//         if (delta < 0 || exp(-delta / temp) > PROB()) {
//             current_inv = NULL;

//             current_inv = copy_invidual(best_neighbour);

//             results[0][*iter_count] = current_inv;

//             (*iter_count)++;

//             tabu_move* best_move = copy_tabu_move(local_tabu_list->moves[best_move_idx]);

//             add_tabu_move(global_tabu_list, best_move);

//             update_tabu_list(global_tabu_list);
//         }

//         delete_tabu_list(local_tabu_list);
//         local_tabu_list = NULL;

//         for (int i = 0; i < neighbour_size; i++)
//             delete_invidual(neighbours[i]);

//         free(neighbours);

//         temp *= alpha;
//     } 

//     delete_hash_table(hs);

//     return results;
// }

invidual*** alg_aoc(problem* prob, int iter_size, int ant_size, int evaporation_flag, double alpha, double beta, double rho) {
    size_t prob_size = prob->distance_matrix.size_x;

    invidual*** results = (invidual***)malloc(iter_size * sizeof(invidual**));

    for (int i = 0; i < iter_size; i++)
        results[i] = (invidual**)malloc(ant_size * sizeof(invidual*));

    for (int i = 0; i < iter_size; i++)
        for (int j = 0; j < ant_size; j++)
            results[i][j] = (invidual*)malloc(sizeof(invidual));

    double** total_pheromones = (double**)malloc(prob_size * sizeof(double*));
    for (int j = 0; j < prob_size; j++)
        total_pheromones[j] = (double*)malloc(prob_size * sizeof(double));

    for (int j = 0; j < prob_size; j++)
        for (int k = 0; k < prob_size; k++)
            total_pheromones[j][k] = 0.0;

    for (int i = 0; i < iter_size; i++) {
        ant** as = (ant**)malloc(ant_size * sizeof(ant*));
        for (int j = 0; j < ant_size; j++) {
            as[j] = (ant*)malloc(sizeof(ant));
            int starting_position = RANDOM(0, prob_size - 1);
            as[j] = create_ant(starting_position, prob_size);
        }

        for (int j = 0; j < ant_size; j++) {
            ant_traverse(as[j], prob->distance_matrix.entries, total_pheromones, alpha, beta, prob_size);
            as[j]->solution->fitness = evaluate_invidual(prob, as[j]->solution);
            update_pheromones(total_pheromones, as[j], prob_size);
            if (evaporation_flag)
                evaporate_pheromones(total_pheromones, rho, prob_size);
        }

        // if (evaporation_flag)
        //     evaporate_pheromones(total_pheromones, rho, prob_size);

        for (int j = 0; j < ant_size; j++)
            results[i][j] = copy_invidual(as[j]->solution);
        
        
        for (int j = 0; j < ant_size; j++) {
            delete_ant(as[j]);
        }
        free(as);
    }

    return results;
}

invidual*** alg_hybrid(problem* prob, int tabu_iter_size, int sa_iter_size, int tabu_neighbour_size, int sa_neightbour_size, int tabu_size, double start_temp, double end_temp, double alpha, int mutation_flag) {
    invidual* current_inv = create_invidual_random(prob->distance_matrix.size_x);
    current_inv->fitness = evaluate_invidual(prob, current_inv);

    invidual*** results = (invidual***)malloc(sizeof(invidual**));

    results[0] = (invidual**)malloc(tabu_iter_size * sizeof(invidual*));

    for (int j = 0; j < tabu_iter_size; j++)
        results[0][j] = (invidual*)malloc(sizeof(invidual));
    
    tabu_list* global_tabu_list = create_tabu_list(tabu_size);
    tabu_list* local_tabu_list = NULL;
    
    void (*m1)(invidual*,int,int);
    void (*m2)(invidual*,int,int);

    if (mutation_flag) {
        m1 = mutation_inversion;
        m2 = mutation_single_swap;
    } else {
        m1 = mutation_single_swap;
        m2 = mutation_inversion;
    }
    
    for (int t = 0; t < tabu_iter_size; t++) {
        if (local_tabu_list == NULL)
            local_tabu_list = create_tabu_list(tabu_size);
        
        
        invidual** neighbours = (invidual**)malloc(tabu_neighbour_size * sizeof(invidual*));
        for (int i = 0; i < tabu_neighbour_size; i++) neighbours[i] = (invidual*)malloc(sizeof(invidual));
        int neighbour_count = 0;
        int best_move_idx = 0;

        tabu_move* temp_holder[tabu_neighbour_size];

        while (neighbour_count < tabu_neighbour_size) {
            int idx1 = RANDOM(0, current_inv->gene.size - 1);
            int idx2 = RANDOM(idx1, current_inv->gene.size - 1);
            tabu_move* tm = create_tabu_move(mutation_flag, idx1, idx2, TABU_TENURE);
            if (!is_tabu(global_tabu_list, tm)) {
                neighbours[neighbour_count] = copy_invidual(current_inv);
                m1(neighbours[neighbour_count], idx1, idx2);
                neighbours[neighbour_count]->fitness = evaluate_invidual(prob, neighbours[neighbour_count]);
                temp_holder[neighbour_count] = tm;
                neighbour_count++;
                // add_tabu_move(global_tabu_list, tm);
            }
        }

        invidual* best_neighbour = NULL;

        for (int j = 0; j < neighbour_count; j++) {
            if (best_neighbour == NULL || neighbours[j]->fitness < best_neighbour->fitness) {
                best_neighbour = neighbours[j];
                best_move_idx = j;
            }
        }

        if (current_inv != NULL) {
            delete_invidual(current_inv);
            current_inv = NULL;
        }

        current_inv = copy_invidual(best_neighbour);

        tabu_move* best_move = copy_tabu_move(temp_holder[best_move_idx]);


        for (int i = 0; i < tabu_neighbour_size; i++)
            delete_invidual(neighbours[i]);

        free(neighbours);

        for (int i = 0; i < tabu_neighbour_size; i++)
            free(temp_holder[i]);

        add_tabu_move(global_tabu_list, best_move);

        double temp = start_temp;
        int iter_count = 0;
        
        while (temp > end_temp && iter_count < sa_iter_size) {
            invidual** neighbours2 = (invidual**)malloc(sa_neightbour_size * sizeof(invidual*));
            for (int i = 0; i < sa_neightbour_size; i++) neighbours2[i] = (invidual*)malloc(sizeof(invidual));
            int neighbour_count2 = 0;
            int best_move_idx2 = 0;

            tabu_move* temp_holder2[sa_neightbour_size];

            while (neighbour_count2 < sa_neightbour_size) {
                int idx1 = RANDOM(0, current_inv->gene.size - 1);
                int idx2 = RANDOM(idx1, current_inv->gene.size - 1);
                tabu_move* tm = create_tabu_move(!mutation_flag, idx1, idx2, TABU_TENURE);
                if (!is_tabu(local_tabu_list, tm)) {
                    neighbours2[neighbour_count2] = copy_invidual(current_inv);
                    m2(neighbours2[neighbour_count2], idx1, idx2);
                    neighbours2[neighbour_count2]->fitness = evaluate_invidual(prob, neighbours2[neighbour_count2]);
                    temp_holder2[neighbour_count2] = tm;
                    neighbour_count2++;
                }
            }

            invidual* best_neighbour2 = NULL;

            for (int j = 0; j < neighbour_count2; j++) {
                if (best_neighbour2 == NULL || neighbours2[j]->fitness < best_neighbour2->fitness) {
                    best_neighbour2 = neighbours2[j];
                    best_move_idx2 = j;
                }
            }

            double delta = best_neighbour2->fitness - current_inv->fitness;
            
            if (delta < 0 || exp(-delta / temp) > PROB()) {
                if (current_inv != NULL) {
                    delete_invidual(current_inv);
                    current_inv = NULL;
                }

                current_inv = copy_invidual(best_neighbour2);

                iter_count++;

                tabu_move* best_move = copy_tabu_move(temp_holder2[best_move_idx2]);

                add_tabu_move(local_tabu_list, best_move);

                update_tabu_list(local_tabu_list);
            }

            for (int i = 0; i < sa_neightbour_size; i++)
                delete_invidual(neighbours2[i]);

            free(neighbours2);

            temp *= alpha;
        }

        if (local_tabu_list != NULL) {
            delete_tabu_list(local_tabu_list);
            local_tabu_list = NULL;
        }

        results[0][t] = copy_invidual(current_inv);
    }

    if (local_tabu_list != NULL)
        delete_tabu_list(local_tabu_list);
    local_tabu_list = NULL;

    delete_invidual(current_inv);
    if (global_tabu_list != NULL)
        delete_tabu_list(global_tabu_list);
    global_tabu_list = NULL;
    
    return results;
}