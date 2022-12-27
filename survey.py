#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import make_dataclass

import spring_counterweight as cw

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
    ("max_error_ratio", float),

    ("arm_length0", float),
    ("spring_origin_x0", float),
    ("spring_origin_y0", float),
    ("cable_length0", float),

    ("arm_length1", float),
    ("spring_origin_x1", float),
    ("spring_origin_y1", float),
    ("cable_length1", float),

    ("spring_od", float),
    ("spring_pn", str)
])
all_solutions = pd.DataFrame(np.zeros(2*len(catalog), dtype=solution_datatype))


# actually do the survey
for spring_qty in range(1,3):
    print(f"Searching the catalog with spring_qty={spring_qty}...")

    for spring_index in range(len(catalog)):
        i = spring_index + (spring_qty-1)*len(catalog)
        item = catalog.iloc[spring_index]
    
        test_spring = cw.Spring(stiffness=item["stiffness"], length_natural=item["length_natural"], length_max=item["length_max"], od=item["od"], pn=item["part_number"])
        soln = cw.spring_optimize(test_spring, spring_qty)
    
        if i % 10 == 0:
            print(f"i = {i} of {len(catalog)}")
    
        # save all the successes to a table
        if soln.success:
            # get the percentage error value
            vc_list = cw.vc_list_from_parameters(soln.x, test_spring, spring_qty)
            net_torque = cw.compute_net_torque(vc_list)
            max_error_ratio = np.max(np.abs(net_torque)) / cw.gravity_torque_max
            
            print("Success", max_error_ratio, soln.x)
            all_solutions.loc[i, "solve_success"] = 1
            all_solutions.loc[i, "spring_index"] = spring_index
            all_solutions.loc[i, "residual"] = soln.fun
            all_solutions.loc[i, "spring_qty"] = spring_qty
            all_solutions.loc[i, "max_error_ratio"] = max_error_ratio

            offset = 4
            all_solutions.loc[i, "arm_length0"] = soln.x[0]
            all_solutions.loc[i, "arm_length1"] = soln.x[0 + offset]

            all_solutions.loc[i, "spring_origin_x0"] = soln.x[1]
            all_solutions.loc[i, "spring_origin_x1"] = soln.x[1 + offset]
        
            all_solutions.loc[i, "spring_origin_y0"] = soln.x[2]
            all_solutions.loc[i, "spring_origin_y1"] = soln.x[2 + offset]
        
            all_solutions.loc[i, "cable_length0"] = soln.x[3]
            all_solutions.loc[i, "cable_length1"] = soln.x[3 + offset]

            all_solutions.loc[i, "spring_od"] = test_spring.od
            all_solutions.loc[i, "spring_pn"] = test_spring.pn

        
# remove all spring solutions that did not successfully solve
all_solutions = all_solutions[all_solutions.solve_success == 1]

# save out the database to an excel file for manual exploration of the solution space
all_solutions.to_csv("solutions.csv")

print("complete")

