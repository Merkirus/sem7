import pandas as pd
import statistics
import matplotlib.pyplot as plt
import re
import numpy as np
from sklearn.preprocessing import MinMaxScaler

GREEDY_VALUES = {
"A-n32-k5": [934.408759,1091.471485, 1064.235604, 1063.815088, 1054.855991, 1069.937979, 1045.571223, 1000.406374,
1132.177379, 956.034743, 1101.284994, 1040.799368, 1085.967667, 1095.288416, 1071.058057, 1034.674520, 992.117317, 1001.865536, 1025.504497, 1003.374971, 1004.090603, 997.588627, 958.806303, 1060.858358,
1071.986475, 1013.720805, 989.515763, 1108.442644, 1028.164372, 1000.566798, 984.156999,1004.993184],
"A-n37-k6": [1166.286367, 1089.416004, 1113.647716, 1121.944824, 1257.323488,
1100.210800, 1247.490040, 1238.532148, 1220.360971, 1211.496743, 1105.265742, 1184.514632, 1065.751265,
1241.317040, 1237.010789, 1174.786296, 1240.374432, 1213.805485, 1081.278475, 1213.090694, 1217.797608, 1119.412564,
1218.434667, 1144.047224, 1223.884238, 1112.124991, 1120.297540, 1106.565801, 1123.705230, 1099.090359, 1168.671924,
1212.881564, 1151.694776, 1117.374015, 1218.983867,
1115.470578, 1099.090359],
"A-n39-k5": [997.832157, 1171.936882, 1077.798285, 1059.482108, 1087.903634, 1036.905829,
1032.662255, 1074.011804, 1016.796597, 1014.530407, 1252.805492, 1018.703607, 986.264013, 1053.303340,
997.832157, 1110.937674, 1062.089595, 1029.377645, 1014.330368, 1135.488018, 1108.180078, 1031.842379, 1054.864261,
1104.192691, 1060.289954, 1114.349227, 1123.352327, 994.174475, 1169.676518, 1120.375261, 984.625812, 1163.541827,
1054.864261, 1145.843358, 1068.476624, 1156.638718,
990.371607, 1081.822867, 1103.776058],
"A-n45-k6": [1168.830835, 1116.802731, 1141.314006, 1118.604280,
1106.739922, 1121.961185, 1175.726383, 1164.752905, 1141.719993, 1303.099460, 1123.861367, 1068.680832, 1244.481325, 1129.943954, 1284.520417, 1209.064374, 1102.177532,
1102.177532, 1102.220324, 1077.133682, 1136.077956, 1182.812045, 1262.900136, 1186.718904, 1213.088335, 1103.814837, 1170.301649, 1169.603260, 1157.655835, 1182.093891, 1070.377342, 1257.495438, 1340.926385, 1338.060623, 1182.567122, 1289.702839,
1239.298782, 1102.909704, 1142.408952, 1177.497798, 1078.367006,
1149.254697, 1099.230127, 1165.294268, 1176.796112],
"A-n48-k7": [1325.218854, 1408.253665,
1395.429475, 1396.302867, 1396.279834, 1465.798671,
1521.097792, 1402.062293, 1442.575477, 1506.714174,
1364.251138, 1404.203971, 1605.831375, 1411.166263, 1440.706116, 1389.358966, 1385.730342, 1440.706116, 1393.188745,
1342.805747, 1394.582128, 1404.963768, 1410.946712, 1421.550995, 1389.973046, 1375.525536, 1392.272069, 1418.484600, 1430.234056, 1393.064421, 1405.085047, 1437.231343, 1440.190956, 1402.897735, 1384.266001, 1366.232945,
1470.395971, 1364.408700, 1337.422469, 1388.717821, 1353.991521, 1499.091975, 1399.077762, 1436.824844,
1394.849187, 1389.362497, 1409.106417, 1388.896428],
"A-n54-k7": [1450.620594, 1624.111690,
1412.601918, 1640.962468,
1615.357794, 1562.951236, 1389.983964, 1431.126599, 1462.686288, 1487.598232, 1491.976996, 1378.818208,
1411.880055, 1379.103332, 1396.810996, 1536.121497, 1418.249393, 1623.581026, 1424.345990, 1487.417281, 1430.930325,
1405.195704, 1638.902604, 1388.235528, 1453.934514, 1606.360759, 1630.061955, 1398.390822, 1555.455139, 1440.479682, 1501.544650, 1468.703077, 1389.717710, 1552.822941,
1466.774386, 1529.794171, 1386.202640, 1475.506100, 1411.243559, 1405.723015, 1410.966527,
1446.020075, 1459.450348, 1458.620594, 1640.962468, 1628.572869, 1497.861974, 1612.045378, 1406.156999,
1502.996947, 1405.377673, 1454.551954,
1501.206553, 1466.743265],
"A-n60-k9": [1762.002729, 1827.523277, 1696.849661, 1749.761367, 1742.360100, 1732.513815, 1661.539336, 1782.481624,
1746.447952, 1783.143460, 1718.552901, 1741.534832, 1721.189073, 1747.762917, 1748.002899, 1755.414510,
1715.353735, 1758.580222, 1828.944467, 1723.148319, 1764.154755, 1742.365856, 1714.954633, 1793.627278,
1763.908626, 1715.353735, 1707.450290, 1754.529754, 1794.040857,
1794.706241, 1728.574854, 1756.295627, 1588.831236,
1813.860446, 1762.628158, 1695.763083, 1598.326836, 1749.530852, 1867.731698, 1754.773029,
1765.288650, 1828.944467, 1757.240722, 1702.569536, 1737.627634, 1706.810336, 1836.182994, 1742.060147,
1787.278571, 1696.911548, 1814.733798, 1599.259855, 1725.731453, 1739.605914,
1748.822123, 1699.380387, 1705.307233,
1748.872461, 1779.731282,
1869.312118]
}

def run(problem, alg, plot_flag):
    results = []

    result_path = f"Solution/{problem}/{problem}-{alg}/results.txt"

    with open(result_path, "w") as f:
        f.write("")

    for i in range(10):
        path = f"Solution/{problem}/{problem}-{alg}/output{i}.csv"

        df = pd.read_csv(path, header=None)

        best_values = []
        worst_values = []
        avg_values = []

        if (alg == "sa" or alg == "hybrid"):
            best_values = df.iloc[0]
        else:
            for index, row in df.iterrows():
                best_values.append(row.min())
                worst_values.append(row.max())
                avg_values.append(row.mean())

        best_value = min(best_values)
        results.append(best_value)
        
        with open(result_path, "a+") as f:
            f.write(f"### Problem: {problem}; Alg: {alg}; Iter: {i}; Best value: {best_value} ###\n")

        if (plot_flag == 1):
            plt.figure(figsize=(10,6))
            
            plt.plot(best_values, label="Best values", linestyle='-', color='green')
            if (alg == "sa" or alg == "hybrid"):
                start_min = 0
                current_min = float('inf')
                size = len(best_values)
                for index, v in enumerate(best_values):
                    if v < current_min:
                        current_min = v
                        temp = index/size
                        plt.axhline(y=v, xmin=start_min, xmax=temp, linestyle='--', color='red')
                        start_min = temp
            else:
                plt.plot(worst_values, label="Worst values", linestyle='-', color='red')
                plt.plot(avg_values, label="Avg values", linestyle='-', color='blue')
            
            plt.xlabel("Iteration")
            plt.ylabel("Value")
            plt.title("Alg plots")
            plt.legend()
            plt.show()
            

    
    with open(result_path, "a+") as f:
        f.write(f"### SUM UP Problem: {problem}; Alg: {alg}; Best value: {min(results)}; Worst value: {max(results)}; Avg value: {sum(results)/len(results)}; Std dev: {statistics.stdev(results)} ###\n")

def compare(problem):
    for i in range(10):
        result_path_ea = f"Solution/{problem}/{problem}-ea/output{i}.csv"
        result_path_ts = f"Solution/{problem}/{problem}-ts/output{i}.csv"
        result_path_sa = f"Solution/{problem}/{problem}-sa/output{i}.csv"
        result_path_aco = f"Solution/{problem}/{problem}-aco/output{i}.csv"
        result_path_hybrid = f"Solution/{problem}/{problem}-hybrid/output{i}.csv"

        best_values_ea = []
        best_values_ts = []
        best_values_sa = []
        best_values_aco = []
        best_values_hybrid = []

        df_ea = pd.read_csv(result_path_ea, header=None)
        df_ts = pd.read_csv(result_path_ts, header=None)
        df_sa = pd.read_csv(result_path_sa, header=None)
        df_aco = pd.read_csv(result_path_aco, header=None)
        df_hybrid = pd.read_csv(result_path_hybrid, header=None)

        best_values_sa = df_sa.iloc[0]
        best_values_hybrid = df_hybrid.iloc[0]

        for index, row in df_ea.iterrows():
            best_values_ea.append(row.min())
        
        for index, row in df_ts.iterrows():
            best_values_ts.append(row.min())

        for index, row in df_aco.iterrows():
            best_values_aco.append(row.min())
        
        plt.figure(figsize=(10,6))
        plt.plot(best_values_ea, label="EA", linestyle='-', color='green')
        plt.plot(best_values_ts, label="TS", linestyle='-', color='blue')
        plt.plot(best_values_sa, label="SA", linestyle='-', color='red')
        plt.plot(best_values_aco, label="ACO", linestyle='-', color='black')
        plt.plot(best_values_hybrid, label="HD", linestyle='-', color='orange')
        plt.xlabel("Iteration")
        plt.ylabel("Value")
        plt.title("Comparison of algorithms")
        plt.legend()
        plt.show()

def compare_normalize(problem):
    for i in range(10):
        result_path_ea = f"Solution/{problem}/{problem}-ea/output{i}.csv"
        result_path_ts = f"Solution/{problem}/{problem}-ts/output{i}.csv"
        result_path_sa = f"Solution/{problem}/{problem}-sa/output{i}.csv"
        result_path_aco = f"Solution/{problem}/{problem}-aco/output{i}.csv"
        result_path_hybrid = f"Solution/{problem}/{problem}-hybrid/output{i}.csv"

        best_values_ea = []
        best_values_ts = []
        best_values_sa = []
        best_values_aco = []
        best_values_hybrid = []

        df_ea = pd.read_csv(result_path_ea, header=None)
        df_ts = pd.read_csv(result_path_ts, header=None)
        df_sa = pd.read_csv(result_path_sa, header=None)
        df_aco = pd.read_csv(result_path_aco, header=None)
        df_hybrid = pd.read_csv(result_path_hybrid, header=None)

        best_values_sa = df_sa.iloc[0]
        best_values_hybrid = df_hybrid.iloc[0]

        for index, row in df_ea.iterrows():
            best_values_ea.append(row.min())
        
        for index, row in df_ts.iterrows():
            best_values_ts.append(row.min())

        for index, row in df_aco.iterrows():
            best_values_aco.append(row.min())
        
        scaler = MinMaxScaler()
        ea_scaled = scaler.fit_transform(np.array(best_values_ea).reshape(-1, 1)).flatten()
        ts_scaled = scaler.fit_transform(np.array(best_values_ts).reshape(-1, 1)).flatten()
        sa_scaled = scaler.fit_transform(np.array(best_values_sa).reshape(-1, 1)).flatten()
        aco_scaled = scaler.fit_transform(np.array(best_values_aco).reshape(-1, 1)).flatten()
        hybrid_scaled = scaler.fit_transform(np.array(best_values_hybrid).reshape(-1, 1)).flatten()

        max_length = max(len(best_values_ea), len(best_values_ts))
        max_length = max(max_length, len(best_values_sa))
        max_length = max(max_length, len(best_values_aco))
        max_length = max(max_length, len(best_values_hybrid))

        ea_int = np.interp(np.linspace(0, len(best_values_ea)-1, max_length), np.arange(len(best_values_ea)), ea_scaled)
        ts_int = np.interp(np.linspace(0, len(best_values_ts)-1, max_length), np.arange(len(best_values_ts)), ts_scaled)
        sa_int = np.interp(np.linspace(0, len(best_values_sa)-1, max_length), np.arange(len(best_values_sa)), sa_scaled)
        aco_int = np.interp(np.linspace(0, len(best_values_aco)-1, max_length), np.arange(len(best_values_aco)), aco_scaled)
        hybrid_int = np.interp(np.linspace(0, len(best_values_hybrid)-1, max_length), np.arange(len(best_values_hybrid)), hybrid_scaled)
        
        plt.figure(figsize=(10,6))
        plt.plot(ea_int, label="EA", linestyle='-', color='green')
        plt.plot(ts_int, label="TS", linestyle='-', color='blue')
        plt.plot(sa_int, label="SA", linestyle='-', color='red')
        plt.plot(aco_int, label="ACO", linestyle='-', color='black')
        plt.plot(hybrid_int, label="HD", linestyle='-', color='orange')
        plt.xlabel("Iteration")
        plt.ylabel("Value")
        plt.title("Comparison of algorithms")
        plt.legend()
        plt.show()

def box(problem):
    result_path_ea = f"Solution/{problem}/{problem}-ea/results.txt"
    result_path_ts = f"Solution/{problem}/{problem}-ts/results.txt"
    result_path_sa = f"Solution/{problem}/{problem}-sa/results.txt"
    result_path_aco = f"Solution/{problem}/{problem}-aco/results.txt"
    result_path_hybrid = f"Solution/{problem}/{problem}-hybrid/results.txt"

    rgx = r"Best value:\s*([0-9.]+)"

    results_ea = []
    results_ts = []
    results_sa = []
    results_aco = []
    results_hybrid = []

    with open(result_path_ea, 'r') as f:
        for i in range(10):
            line = f.readline()
            if line:
                mtch = re.search(rgx, line)
                if mtch:
                    results_ea.append(float(mtch.group(1)))

    with open(result_path_ts, 'r') as f:
        for i in range(10):
            line = f.readline()
            if line:
                mtch = re.search(rgx, line)
                if mtch:
                    results_ts.append(float(mtch.group(1)))

    with open(result_path_sa, 'r') as f:
        for i in range(10):
            line = f.readline()
            if line:
                mtch = re.search(rgx, line)
                if mtch:
                    results_sa.append(float(mtch.group(1)))
    
    with open(result_path_aco, 'r') as f:
        for i in range(10):
            line = f.readline()
            if line:
                mtch = re.search(rgx, line)
                if mtch:
                    results_aco.append(float(mtch.group(1)))

    with open(result_path_hybrid, 'r') as f:
        for i in range(10):
            line = f.readline()
            if line:
                mtch = re.search(rgx, line)
                if mtch:
                    results_hybrid.append(float(mtch.group(1)))

    data = [
        GREEDY_VALUES[problem],
        results_ea,
        results_ts,
        results_sa,
        results_aco,
        results_hybrid
    ]

    plt.figure(figsize=(8,6))
    box = plt.boxplot(data, patch_artist=True, vert=True)
    plt.title("Comparison of algorithms")
    colors = ["yellow", "green", "red", "blue", "black", "orange"]
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    plt.xticks([1,2,3,4,5,6], ["GR", "EA", "TS", "SA", "ACO", "HD"])
    plt.ylabel("Value")

    plt.show()

def best_route(problem):
    best_route = None
    best_route_value = float('inf')

    for i in range(10):
        result_path_sa = f"Solution/{problem}/{problem}-hybrid/output{i}.csv"

        df_sa = pd.read_csv(result_path_sa, header=None)

        best_values_sa = df_sa.iloc[0]

        min_value = float('inf')
        min_index = 0

        for index, value in enumerate(best_values_sa):
            if value < min_value:
                min_value = value
                min_index = index

        result_path_sa = f"Solution/{problem}/{problem}-hybrid/str_output{i}.csv"
        df_sa_str = pd.read_csv(result_path_sa, header=None)

        genes = df_sa_str.iloc[0]

        if min_value < best_route_value:
            best_route = genes[min_index]
        
    print(best_route)