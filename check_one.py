#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import spring_counterweight as cw

catalog = pd.read_csv("springs_english.csv", sep=",")


########
# Choose a specific spring solution of interest. This should already be optimized.
spring_qty = 1
index = 1051 # index in the spring catalog
soln = np.array([ 0.0176,  0.0427, -0.2117,  0.000,
        1.50000000e-02, -8.60161725e-02, -2.70000000e-01,  3.28437002e-17])

# single spring best
# spring_qty = 1
# index = 921
# soln = np.array([0.015, 0.024, -0.265, 0.000, 
#                 0.015, 0.050, -0.200, 0.000])

## Now actually run the calculation
item = catalog.iloc[index]
test_spring = cw.Spring(stiffness=item["stiffness"], length_natural=item["length_natural"], length_max=item["length_max"], od=item["od"], pn=item["part_number"], force_init=item["load_min"])
print(test_spring)


cw.show_solution(test_spring, spring_qty, soln)

