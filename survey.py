#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import make_dataclass

import simple_spring as cw

print("Loading the spring catalog...")
catalog = pd.read_csv("springs_english.csv", sep=",")
print(f"The spring catalog contains {len(catalog)} options")
print(catalog)

# setup the list of things to record for each solution
solution_datatype = np.dtype([
    ("spring_index", int),
    ("residual", float),
    ("spring_qty", int),
    ("solve_success", int),

    ("arm_length0", float),
    ("spring_origin_x0", float),
    ("spring_origin_y0", float),
    ("cable_length0", float),

    ("arm_length1", float),
    ("spring_origin_x1", float),
    ("spring_origin_y1", float),
    ("cable_length1", float),
])

all_solutions = pd.DataFrame(np.zeros(len(catalog), dtype=solution_datatype))


best_value = 1e9
best_index = -1
best_soln = None

# TODO go back to searching the entire thing
for i in range(240, len(catalog)):
    item = catalog.iloc[i]
    
    spring_qty = 1
    test_spring = cw.Spring(stiffness=item["stiffness"], length_natural=item["length_natural"], length_max=item["length_max"])
    soln = cw.spring_optimize(test_spring, spring_qty)
    
    print("i =", i)

    # track the very best one
    if soln.success and soln.fun < best_value:
        print("New best.")#, i, soln)
        best_value = soln.fun
        best_index = i
        best_soln = soln
    
    # save all the successes to a table
    if soln.success:
        print(soln)
        all_solutions.loc[i, "solve_success"] = 1
        all_solutions.loc[i, "spring_index"] = i
        all_solutions.loc[i, "residual"] = soln.fun
        all_solutions.loc[i, "spring_qty"] = spring_qty

        offset = 4
        all_solutions.loc[i, "arm_length0"] = soln.x[0]
        all_solutions.loc[i, "arm_length1"] = soln.x[0 + offset]

        all_solutions.loc[i, "spring_origin_x0"] = soln.x[1]
        all_solutions.loc[i, "spring_origin_x1"] = soln.x[1 + offset]
        
        all_solutions.loc[i, "spring_origin_y0"] = soln.x[2]
        all_solutions.loc[i, "spring_origin_y1"] = soln.x[2 + offset]
        
        all_solutions.loc[i, "cable_length0"] = soln.x[3]
        all_solutions.loc[i, "cable_length1"] = soln.x[3 + offset]

        
# remove all spring solutions that did not successfully solve
all_solutions = all_solutions[all_solutions.solve_success == 1]

# TODO save out the database to an excel file for manual exploration of the solution space
all_solutions.to_csv("solutions.csv")

print("complete")

