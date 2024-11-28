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


PROBLEM1 = "A-n32-k5"
PROBLEM2 = "A-n37-k6"
PROBLEM3 = "A-n39-k5"
PROBLEM4 = "A-n45-k6"
PROBLEM5 = "A-n48-k7"
PROBLEM6 = "A-n54-k7"
PROBLEM7 = "A-n60-k9"

PROBLEMS = [
    PROBLEM1,
    PROBLEM2,
    PROBLEM3,
    PROBLEM4,
    PROBLEM5,
    PROBLEM6,
    PROBLEM7
]

ITER_SIZE = 10

# gen_size, pop_size, elit_flag, select_flag, cross_flag, mut_flag, px, pm
EA_ARGS = [10000, 50, 0, 2, 1, 0, 0.2, 0.4]
# iter_size, neigh_size, tabu_size, mut_flag, greedy_flag
TS_ARGS = [10000, 20, 100, 0, 0]
# iter_size, neigh_size, tabu_size, temp_start, temp_end, alpha, mut_flag
SA_ARGS = [50000, 3, 100, 100.0, 0.001, 0.99995, 0]

def move_files(alg, problem):
    dest_path = f"Solution/{problem}/{problem}-{alg}"

    if not os.path.exists(dest_path):
        os.makedirs(dest_path)
    
    for i in range(10):
        file = f"output{i}.csv"
        dest_file = os.path.join(dest_path, file)
        shutil.move(file, dest_file)
        file = f"str_output{i}.csv"
        dest_file = os.path.join(dest_path, file)
        shutil.move(file, dest_file)

for p in PROBLEMS:
    p_path = f"Prob/{p}.vrp".encode('utf-8')
    # ts1 = time.perf_counter()
    # metaheuristics.run_ea(p_path, EA_ARGS[0], EA_ARGS[1], EA_ARGS[2], EA_ARGS[3], EA_ARGS[4], EA_ARGS[5], EA_ARGS[6], EA_ARGS[7])
    # metaheuristics.run_ts(p_path, TS_ARGS[0], TS_ARGS[1], TS_ARGS[2], TS_ARGS[3], TS_ARGS[4])
    # metaheuristics.run_sa(p_path, SA_ARGS[0], SA_ARGS[1], SA_ARGS[2], SA_ARGS[3], SA_ARGS[4], SA_ARGS[5], SA_ARGS[6])
    # te1 = time.perf_counter()
    # move_files("ea", p)
    # move_files("ts", p)
    # move_files("sa", p)
    # run(p, "ea", 1)
    # run(p, "ts", 0)
    # run(p, "sa", 0)
    # print(f"Exec time: {te1 - ts1} sec.")
    # box(p)
    # compare(p)
    # compare_normalize(p)
    # best_route(p)
