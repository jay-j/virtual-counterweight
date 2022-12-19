#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import simple_spring as cw

catalog = pd.read_csv("springs_english.csv", sep=",")

spring_qty = 2

# Choose a specific spring solution of interest. This should already be optimized.
index = 1010 # index in the spring catalog
soln = np.array([ 1.51075338e-02,  1.30000000e-01, -2.70000000e-01,  0.00000000e+00,
        1.50000000e-02, -8.60161725e-02, -2.70000000e-01,  3.28437002e-17])

# double spring exploring
index = 952
soln = np.array([0.020, 0.130, -0.270, 0.000, 
                0.019, -0.082, -0.270, 0.000])

# single spring best
# index = 835
# soln = np.array([0.021, 0.032, -0.231, 0.000, 
#                 0.015, 0.050, -0.200, 0.000])
# spring_qty = 1

item = catalog.iloc[index]
test_spring = cw.Spring(stiffness=item["stiffness"], length_natural=item["length_natural"], length_max=item["length_max"])

cw.show_solution(test_spring, spring_qty, soln)

