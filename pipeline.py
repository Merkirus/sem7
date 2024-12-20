import ctypes
import shutil
import os
import time
from plot import run, box, compare, compare_normalize, best_route

metaheuristics = ctypes.CDLL('./lib/metaheuristics.so')
#                                   filename        gen_size        pop_size    elit_flag       select_flag cross_flag      mut_flag        px              pm
metaheuristics.run_ea.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double]
metaheuristics.run_ea.restype = None

#                                   filename        iter_size      neigh_size    tabu_size     mut_flag     greedy_flag
metaheuristics.run_ts.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int]
metaheuristics.run_ts.restype = None

#                                   filename        iter_size      neigh_size    tabu_size     temp_start     temp_end          alpha              mut_flag
metaheuristics.run_sa.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int]
metaheuristics.run_sa.restype = None

#                                   filename        iter_size      ant_size    evapor_flag     alpha           beta               rho
metaheuristics.run_aco.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double]
metaheuristics.run_aco.restype = None

#                                   filename        tabu_iter_size      sa_iter    tabu_neig     sa_neigh           tabu_size   start_temp       end_temp      alpga           mut
metaheuristics.run_hybrid.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int]
metaheuristics.run_hybrid.restype = None

PROBLEM1 = "A-n32-k5"
PROBLEM2 = "A-n37-k6"
PROBLEM3 = "A-n39-k5"
PROBLEM4 = "A-n45-k6"
PROBLEM5 = "A-n48-k7"
PROBLEM6 = "A-n54-k7"
PROBLEM7 = "A-n60-k9"

PROBLEMS = [
    # PROBLEM1,
    # PROBLEM2,
    # PROBLEM3,
    PROBLEM4,
    # PROBLEM5,
    # PROBLEM6,
    # PROBLEM7
]

ITER_SIZE = 10

# gen_size, pop_size, elit_flag, select_flag, cross_flag, mut_flag, px, pm
EA_ARGS = [500, 10000, 50, 5, 1, 0, 0.9, 0.1]
# iter_size, neigh_size, tabu_size, mut_flag, greedy_flag
TS_ARGS = [10000, 20, 100, 0, 0]
# iter_size, neigh_size, tabu_size, temp_start, temp_end, alpha, mut_flag
SA_ARGS = [50000, 3, 100, 100.0, 0.001, 0.99995, 0]
# iter, ants, evapor, alha, beta, rho
ACO_ARGS = [200, 50, 1, 0.55, 5.0, 0.10]
#tabu_iter_size, sa_iter, tabu_neig, sa_neigh, tabu_size, start_temp, end_temp, alpha           mut
HYBRID_ARGS = [500, 5000, 20, 5, 100, 100.0, 0.001, 0.999, 1]

def move_files(alg, problem):
    dest_path = f"Solution/{problem}/{problem}-{alg}"

    if not os.path.exists(dest_path):
        os.makedirs(dest_path)
    
    for i in range(ITER_SIZE):
        file = f"output{i}.csv"
        dest_file = os.path.join(dest_path, file)
        shutil.move(file, dest_file)
        file = f"str_output{i}.csv"
        dest_file = os.path.join(dest_path, file)
        shutil.move(file, dest_file)

for p in PROBLEMS:
    p_path = f"Prob/{p}.vrp".encode('utf-8')
    # ts1 = time.perf_counter()
    # metaheuristics.run_aco(p_path, ACO_ARGS[0], ACO_ARGS[1], ACO_ARGS[2], ACO_ARGS[3], ACO_ARGS[4], ACO_ARGS[5])
    # move_files("aco", p)
    # run(p, "aco", 1)
    # metaheuristics.run_hybrid(p_path, HYBRID_ARGS[0], HYBRID_ARGS[1], HYBRID_ARGS[2], HYBRID_ARGS[3], HYBRID_ARGS[4], HYBRID_ARGS[5], HYBRID_ARGS[6], HYBRID_ARGS[7], HYBRID_ARGS[8])
    # move_files("hybrid", p)
    # run(p, "hybrid", 1)
    # metaheuristics.run_ea(p_path, EA_ARGS[0], EA_ARGS[1], EA_ARGS[2], EA_ARGS[3], EA_ARGS[4], EA_ARGS[5], EA_ARGS[6], EA_ARGS[7])
    # metaheuristics.run_ts(p_path, TS_ARGS[0], TS_ARGS[1], TS_ARGS[2], TS_ARGS[3], TS_ARGS[4])
    # metaheuristics.run_sa(p_path, SA_ARGS[0], SA_ARGS[1], SA_ARGS[2], SA_ARGS[3], SA_ARGS[4], SA_ARGS[5], SA_ARGS[6])
    # te1 = time.perf_counter()
    # move_files("ea", p)
    # move_files("ts", p)
    # move_files("sa", p)
    # run(p, "ea", 0)
    # run(p, "ts", 0)
    # run(p, "sa", 0)
    # print(f"Exec time: {te1 - ts1} sec.")
    # box(p)
    # compare(p)
    # compare_normalize(p)
    best_route(p)
