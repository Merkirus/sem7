#include <string>
#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <random>
#include <time.h>

#define RANDOM(min, max) (rand() % ((max) - (min) + 1) + (min))

using namespace std;

double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow((x2-x1), 2) + pow((y2-y1), 2));
}

class Problem {
    public:
        string name;
        string comment;
        string type;
        int dim;
        string edge_weight_type;
        int cap;
        double** d_matrix;
        int* demand;
        Problem(string name, string comment, string type, int dim, string edge_weight_type, int cap) {
            this->name = name;
            this->comment = comment;
            this->type = type;
            this->dim = dim;
            this->edge_weight_type = edge_weight_type;
            this->cap = cap;
        }
        ~Problem() {
            for (int i = 0; i < this->dim; i++) {
                delete[] d_matrix[i];
            }
            delete[] d_matrix;
            delete[] demand;
        }
        void set_d_matrix(double** d_matrix) {
            this->d_matrix = d_matrix;
        }
        void set_demand(int* demand) {
            this->demand = demand;
        }
        void print_d_matrix() {
            for (int i = 0; i < this->dim; i++) {
                for (int j = 0; j < this->dim; j++) {
                    cout << this->d_matrix[i][j] << " ";
                }
                cout << endl;
            }
        }
        void print_demand() {
            for (int i = 0; i < this->dim; i++) {
                cout << this->demand[i] << " ";
            }
            cout << endl;
        }
};

class Loader {
    public:
        Loader(string file_path) {
            this->file_path = file_path;
            this->get_extension_from_file();

            if (this->extension == "vrp") {
                this->prob = this->load_vrb_problem();
            }
        }
        string get_extension() {
            return this->extension;
        }
        Problem* get_problem() {
            return prob;
        }
    private:
        const int META_VRP_SIZE = 6;
        string file_path;
        string extension;
        Problem* prob;
        void get_extension_from_file() {
            regex reg(R"(\.(.*))");
            smatch matches;

            if (regex_search(file_path, matches, reg)) {
                if (matches.size() > 1) {
                    this->extension = matches[1];
                }
            }
        }
        Problem* load_vrb_problem() {
            ifstream file(this->file_path);

            if (file.is_open()) {

                string line;

                vector<string> vec_meta;

                for (int i = 0; i < META_VRP_SIZE; i++) {
                    getline(file, line);
                    vec_meta.push_back(this->get_metadata(line));
                }

                Problem* prob = new Problem(vec_meta.front(),vec_meta.at(1),vec_meta.at(2),stoi(vec_meta.at(3)),vec_meta.at(4),stoi(vec_meta.back()));

                int coords[prob->dim][2];

                getline(file, line); // COORDS
                
                for (int i = 0; i < prob->dim; i++) {
                    getline(file, line);
                    istringstream iss(line);
                    vector<int> vec_coords;
                    
                    int value;
                    while (iss >> value) vec_coords.push_back(value);
                    coords[i][0] = vec_coords.at(1);
                    coords[i][1] = vec_coords.at(2);
                }

                getline(file, line); // DEMAND_SECTION

                int* demand = new int[prob->dim];

                for (int i = 0; i < prob->dim; i++) {
                    getline(file, line);
                    istringstream iss(line);
                    vector<int> vec_demand;

                    int value;
                    while (iss >> value) vec_demand.push_back(value);
                    demand[i] = vec_demand.at(1);
                }

                double** d_matrix = new double*[prob->dim];
                for (int i = 0; i < prob->dim; i++) d_matrix[i] = new double[prob->dim];

                for (int i = 0; i < prob->dim; i++) {
                    for (int j = 0; j < prob->dim; j++) {
                        int x1 = coords[i][0];
                        int y1 = coords[i][1];
                        int x2 = coords[j][0];
                        int y2 = coords[j][1];
                        d_matrix[i][j] = distance(x1,y1,x2,y2);
                    }
                }
                prob->set_d_matrix(d_matrix);
                prob->set_demand(demand);

                return prob;
            }
            return nullptr;
        }
        string get_metadata(string& line) {
            regex regex(R"(\w+\s*:\s*(.*))");
            smatch matches;

            if (regex_search(line, matches, regex)) {
                if (matches.size() > 1) {
                    return matches[1];
                }
            }

            return "";
        }
};

class Invidual {
    public:
    Invidual() {}
    Invidual(vector<int>&& v) {
        this->result = move(v);
    }
        vector<int> result;
        void gen_random(int n, int depot_number) {
            random_device rd;
            mt19937 eng(rd());

            for (int i = 1; i < n; i++) result.push_back(i);

            shuffle(result.begin(), result.end(), eng);
        }
};

class Method {
    public:
        Method(Problem* prob, Invidual* inv) {
            
        }
};

double test(Problem* prob, Invidual* inv) {    
    int curr_cap = prob->cap;
    int curr_station = 0; // depot
    double d = 0.0;
    double d_depot = 0.0;
    vector<int> route;

    for (int i = 0; i < inv->result.size(); i++) {
        int next_station = inv->result.at(i);
        int curr_demand = prob->demand[next_station];
        if (curr_cap - curr_demand < 0) {
            d_depot += prob->d_matrix[curr_station][0];
            curr_cap = prob->cap;
            curr_station = 0;
            i--;
            route.push_back(0);
            continue;
        } 
        curr_cap -= curr_demand;
        d += prob->d_matrix[curr_station][next_station];
        route.push_back(next_station);
        curr_station = next_station;
    }

    // cout << "ROUTE: ";

    // for (auto& r : route) {
    //     cout << r << " ";
    // }


    // cout << endl;

    d_depot += prob->d_matrix[curr_station][0];

    return d + d_depot;
}
int exists_in_subarray(int* array, int start, int end, int value) {
    for (int i = start; i <= end; i++) {
        if (array[i] == value) {
            return 1; // Found
        }
    }
    return 0; // Not found
}
Invidual* crossover_ox(Invidual* inv1, Invidual* inv2) {
    size_t inv_size = inv1->result.size();

    int* child = new int[inv_size];

    int start_sub = RANDOM(0, inv_size-1);
    int stop_sub = RANDOM(start_sub, inv_size-1);
    
    int sub_size = stop_sub - start_sub + 1;

    for (int i = 0; i < inv_size; i++) {
        child[i] = -1; // Initialize the child with invalid values
    }
    for (int i = start_sub; i <= stop_sub; i++) {
        child[i] = inv1->result[i];
    }

    // Step 3: Fill the rest of the child with parent2's elements, maintaining the order
    int current = (stop_sub + 1) % inv_size;  // Start filling from the end of the second crossover point
    for (int i = 0; i < inv_size; i++) {
        int parent2_value = inv2->result[(stop_sub + 1 + i) % inv_size];
        if (!exists_in_subarray(child, start_sub, stop_sub, parent2_value)) {
            child[current] = parent2_value;
            current = (current + 1) % inv_size; // Move to the next position, wrapping around
        }
    }

    return new Invidual(vector<int>(child, child + inv_size));
}

void mutation_swap(Invidual* inv) {
    for (int i = 0; i < inv->result.size(); i++) {
        for (int j = i + 1; j < inv->result.size(); j++) {
            if (rand() < 0.1) {
                int temp = inv->result[i];
                inv->result[i] = inv->result[j];
                inv->result[j] = temp;
            }
        }
    }
}

int main(int argc, char const *argv[])
{
    srand(time(NULL));
    Loader loader("Prob/test.vrp");
    // Invidual spec{};
    // Invidual spec2{};
    // Invidual* spec3 = crossover_ox(&spec,&spec2);
    Problem* prob;
    prob = loader.get_problem();

    int pop_size = 100;
    int iters = 1000;

    Invidual** pops = new Invidual*[iters];
    for (int i = 0; i < iters; i++) {
        pops[i] = new Invidual[pop_size];
    }


    for (int i = 0; i < pop_size; i++) {
        pops[0][i].gen_random(prob->dim, 0);
    }

    double results[iters];
    for (int i = 0; i < iters; i++) {
        results[i] = 0.0;
    }

    for (int pop_iter = 0; pop_iter < iters - 1; pop_iter++) {
        for (int j = 0; j < pop_size; j++) {
            int rand_idx = RANDOM(0,pop_size-1);
            int rand_idx2 = RANDOM(0,pop_size-1);
            Invidual* curr_pop = pops[pop_iter];
            Invidual inv1 = curr_pop[rand_idx];
            Invidual inv2 = curr_pop[rand_idx2];
            Invidual* os;
            if (rand() < 0.7) {
                os = crossover_ox(&inv1, &inv2);
            } else {
                os = &inv1;
            } 
            mutation_swap(os);
            pops[pop_iter+1][j].result = os->result;
            double result = test(prob, os);
            cout << result << endl;
            if (results[pop_iter+1] > result) {
                results[pop_iter+1] = result;
            }
        }
    }
    // prob->print_d_matrix();
    // prob->print_demand();
    // spec.gen_random(prob->dim, 0);
    // for (auto& i : spec.result) {
    //     cout << i << " ";
    // }
    // cout << endl;
    // cout << test(prob, &spec) << endl;
    // if (prob != nullptr) {
    //     delete prob;
    // }
    // for (auto& i : results) {
    //     cout << i << endl;
    // }
    return 0;
}
