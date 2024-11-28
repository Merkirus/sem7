#include "structure.h"

#define RANDOM(min, max) (rand() % ((max) - (min) - (min) + 1) + (min))
#define BLANK -100
#define BUFFER 256
#define METADATA_SIZE 7
#define DEPOT 0

invidual* create_invidual_empty(size_t gene_size) {
    invidual* inv = (invidual*)malloc(sizeof(invidual));
    inv->gene.size = gene_size;
    inv->gene.entries = (int*)malloc(gene_size * sizeof(int));
    return inv;
}

invidual* create_invidual_random(size_t gene_size) {
    invidual* inv = (invidual*)malloc(sizeof(invidual));
    inv->gene.size = gene_size;
    inv->gene.entries = (int*)malloc(gene_size * sizeof(int));
    
    for (int i = 0; i < gene_size; i++)
        inv->gene.entries[i] = i;
    
    for (int i = gene_size - 1; i > 0; i--) {
        int j = rand() % (i + 1);

        int temp = inv->gene.entries[i];
        inv->gene.entries[i] = inv->gene.entries[j];
        inv->gene.entries[j] = temp;
    }

    return inv;
}

invidual* create_invidual_greedy(int start, int** closest_stations, size_t gene_size) {
    invidual* inv = (invidual*)malloc(sizeof(invidual));
    inv->gene.size = gene_size;
    inv->gene.entries = (int*)malloc(gene_size * sizeof(int));

    inv->gene.entries[0] = start;
    int current = start;

    int visited_stations[gene_size];
    for (int i = 0; i < gene_size; i++) visited_stations[i] = 0; // false
    visited_stations[current] = 1;

    for (int i = 1; i < gene_size; i++) {
        for (int j = 0; j < gene_size; j++) {
            int next_station = closest_stations[current][j];
            if (visited_stations[next_station]) continue;
            visited_stations[next_station] = 1;
            current = next_station;
            break;
        }
        inv->gene.entries[i] = current;
    }

    return inv;
}

invidual* copy_invidual(invidual* inv) {
    invidual* new_inv = (invidual*)malloc(sizeof(invidual));
    new_inv->fitness = inv->fitness;

    new_inv->gene.size = inv->gene.size;
    new_inv->gene.entries = (int*)malloc(new_inv->gene.size * sizeof(int));

    for (int i = 0; i < inv->gene.size; i++)
        new_inv->gene.entries[i] = inv->gene.entries[i];

    return new_inv;
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
    int curr_station = DEPOT;

    double distance = 0.0;
    int route[BUFFER];
    for (int i = 0; i < BUFFER; i++) route[i] = BLANK;

    for (int i = 0, route_index = 0; i < inv->gene.size; i++) {
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

tabu_move* create_tabu_move(int move_type, int idx1, int idx2, int lifespan) {
    tabu_move* tm = (tabu_move*)malloc(sizeof(tabu_move));
    tm->idx1 = idx1;
    tm->idx2 = idx2;
    tm->move_type = move_type;
    tm->lifespan = lifespan;
    return tm;
}

tabu_move* copy_tabu_move(tabu_move* tm) {
    tabu_move* new_tm = (tabu_move*)malloc(sizeof(tabu_move));
    new_tm->idx1 = tm->idx1;
    new_tm->idx2 = tm->idx2;
    new_tm->move_type = tm->move_type;
    new_tm->lifespan = tm->lifespan;
    return new_tm;
}

tabu_list* create_tabu_list(size_t tabu_size) {
    tabu_list* tl = (tabu_list*)malloc(sizeof(tabu_list));

    tl->moves_size = tabu_size;
    tl->moves = (tabu_move**)malloc(tabu_size * sizeof(tabu_move*));
    for (int i = 0; i < tabu_size; i++)
        tl->moves[i] = (tabu_move*)malloc(sizeof(tabu_move));

    for (int i = 0; i < tabu_size; i++)
        tl->moves[i] = NULL;

    tl->current_size = 0;

    return tl;
}

void add_tabu_move(tabu_list* tl, tabu_move* tm) {
    if (tl->current_size == tl->moves_size) {
        return;
    }
    tl->moves[tl->current_size] = tm;
    tl->current_size++;
}

int is_tabu(tabu_list* tl, tabu_move* tm) {
    for (int i = 0; i < tl->current_size; i++) {
        tabu_move* curr_move = tl->moves[i];
        if (curr_move->idx1 == tm->idx1 &&
            curr_move->idx2 == tm->idx2 &&
            curr_move->move_type == tm->move_type) {
                return 1;
            }
    }
    return 0;
}

void update_tabu_list(tabu_list* tl) {
    for (int i = 0; i < tl->current_size; i++) {
        tl->moves[i]->lifespan--;
    }
    for (int i = 0; i < tl->current_size; i++) {
        if (tl->moves[i]->lifespan <= 0) {
            for (int j = i; j < tl->current_size - 1; j++) {
                tl->moves[j] = tl->moves[j+1];
            }
            tl->current_size--;
            i--;
        }
    }
}

void delete_tabu_list(tabu_list* tl) {
    if (tl != NULL) {
        if (tl->moves != NULL) {
            for (int i = 0; i < tl->current_size; i++)
                if (tl->moves[i] != NULL) {
                    free(tl->moves[i]);
                    tl->moves[i] = NULL;
                }
            free(tl->moves);
            tl->moves = NULL;
        }
        free(tl);
    }
}

problem* create_problem(char* file_path) {
    problem* prob = (problem*)malloc(sizeof(problem));
    FILE* file;
    char line[BUFFER];

    file = fopen(file_path, "r");

    if (file == NULL)
        return NULL;

    for (int i = 0; i < METADATA_SIZE; i++) {
        if (!fgets(line, sizeof(line), file))
            return NULL;

        // char* value = _get_file_data(line);
        regex_t regex;
        regmatch_t matches[2];
        char* pattern = "\\w+\\s*:\\s*(.*)";
        char* value = NULL;
        
        if (regcomp(&regex, pattern, REG_EXTENDED) != 0)
            return NULL;

        if (regexec(&regex, line, 2, matches, 0) == 0) {
            if (matches[1].rm_so != -1) {
                int length = matches[1].rm_eo - matches[1].rm_so;
                value = (char*)malloc(length + 1);
                if (value) {
                    strncpy(value, line + matches[1].rm_so, length);
                    value[length] = '\0';
                }
            }
        }
        
        regfree(&regex);
        
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

    prob->distance_matrix.entries = (double**)malloc(prob->distance_matrix.size_x * sizeof(double*));
    for (int i = 0; i < prob->distance_matrix.size_x; i++)
        prob->distance_matrix.entries[i] = (double*)malloc(prob->distance_matrix.size_y * sizeof(double));

    prob->demand.entries = (int*)malloc(prob->demand.size * sizeof(int));
    
    int coordinates[prob->distance_matrix.size_x][2];

    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        if (!fgets(line, sizeof(line), file))
            return NULL;

        if (!sscanf(line, "%*d %d %d", &coordinates[i][0], &coordinates[i][1]))
            return NULL;
    }

    if (!fgets(line, sizeof(line), file)) // Demand section
        return NULL;

    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        if (!fgets(line, sizeof(line), file))
            return NULL;
        
        if (!sscanf(line, "%*d %d", &prob->demand.entries[i]))
            return NULL;
    }

    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        for (int j = 0; j < prob->distance_matrix.size_y; j++) {
            int x1 = coordinates[i][0];
            int y1 = coordinates[i][1];
            int x2 = coordinates[j][0];
            int y2 = coordinates[j][1];
            double distance = sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
            prob->distance_matrix.entries[i][j] = distance;
        }
    }

    fclose(file);

    return prob;
}

void delete_problem(problem* prob) {
    free(prob->demand.entries);
    
    for (int i = 0; i < prob->distance_matrix.size_x; i++) {
        free(prob->distance_matrix.entries[i]);
    }
    free(prob->distance_matrix.entries);
}