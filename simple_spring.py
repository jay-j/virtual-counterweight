#!/usr/bin/env python3
import numpy as np
import scipy.optimize as optim
import matplotlib.pyplot as plt

## DESIGN PARAMETERS #########################################################

g = 9.81 # acceleration of gravity, m/s^2
mass = 0.30 # kg, mass of the end effector
L = 0.300   # m, length from pivot to the 

# angle of the main arm, radians
# measured with positive DOWN, straight horizontal zero
theta = np.linspace(-np.radians(70), np.radians(70), 32)
theta_label = "Arm Angle, theta (radians)"

# uncompensated torque to support the arm, Nm
gravity_torque = mass*g * L * np.cos(theta)
gravity_torque_max = np.max(gravity_torque)

plt.figure(1)
plt.plot(theta, gravity_torque, label="Uncompensated gravity")
plt.xlabel(theta_label)
plt.ylabel("Torque, Nm")
plt.grid(True)

def rotate_z(vector, angle_rad):
    # one vector, many angles
    vector = vector.reshape((-1,3))
    result = np.zeros((angle_rad.shape[0], 3))
    result[:,0] = vector[:,0]*np.cos(angle_rad) - vector[:,1]*np.sin(angle_rad)
    result[:,1] = vector[:,0]*np.sin(angle_rad) + vector[:,1]*np.cos(angle_rad)
    result[:,2] = vector[:,2]
    return result


# use a model where the lever parametrization is along +X only. Can rotate the entire system later if needed
class VirtualCounterweight:
    def __init__(self, lever_radius, spring_origin, spring_cable_length, spring_stiffness, spring_length_natural, spring_length_max):
        self.lever = np.array([lever_radius, 0, 0])
        self.spring_origin = spring_origin
        self.spring_stiffness = spring_stiffness
        self.spring_length_natural = spring_length_natural
        self.spring_length_max = spring_length_max
        self.static_length = spring_cable_length # 0.010
    
    def compute_spring_length(self, theta):
        lever = rotate_z(self.lever, theta)
        spring_vec = self.spring_origin - lever
        spring_length = np.linalg.norm(spring_vec, axis=1).reshape((-1,1))
        return spring_length
    

    def compute_torque(self, theta):
        # compute new endpoint given radius
        lever = rotate_z(self.lever, theta)

        spring_vec = self.spring_origin - lever
        
        # compute spring direction
        spring_dir = spring_vec / np.linalg.norm(spring_vec, axis=1).reshape((-1,1))

        # compute spring length
        spring_length_delta = np.linalg.norm(spring_vec, axis=1).reshape((-1,1)) - self.static_length - self.spring_length_natural
        # print("spring length delta", spring_length_delta)

        # compute spring force
        spring_force = spring_length_delta * self.spring_stiffness * spring_dir
        # print("spring force:", spring_force)

        # compute torque = r cross F
        torque = lever[:,0]*spring_force[:,1] + lever[:,1]*spring_force[:,0]
        return torque
    



# TODO choose between different real spring options
# should unique spring stiffnesses per spring be allowed? 
spring_stiffness = 450
spring_qty = 2
spring_length_natural = 0.100
spring_length_max = 0.200
parms_per_spring = 4

def compute_net_torque(vc_list):
    net_torque = mass*g * L * np.cos(theta)
    for vc in vc_list:
        net_torque += vc.compute_torque(theta)
    return net_torque


def residual(parm):
    # the function to optimize!
    # six parameters, dimensional locations of the spring
    # given a constant spring

    vc_list = []
    for i in range(spring_qty):
        offset = parms_per_spring * i
        vc_list.append(VirtualCounterweight(lever_radius=parm[offset+0],
            spring_origin=np.array([parm[offset+1], parm[offset+2], 0]),
            spring_cable_length=parm[offset+3],
            spring_stiffness=spring_stiffness,
            spring_length_natural=spring_length_natural,
            spring_length_max=spring_length_max))
    tau = compute_net_torque(vc_list)

    # use least squares error
    error = np.sum(np.power(tau, 2))
    return error
    
def constraint_spring_minimum(parm):
    # inequality type constraint. Valid when return > 0
    # spring delta length must be > 0 (always in tension)
    minimum_margin = 1
    for i in range(spring_qty):
        offset = parms_per_spring * i
        vc = VirtualCounterweight(lever_radius=parm[offset+0],
                spring_origin=np.array([parm[offset+1], parm[offset+2], 0]),
                spring_cable_length=parm[offset+3],
                spring_stiffness=spring_stiffness,
                spring_length_natural=spring_length_natural,
                spring_length_max=spring_length_max)
        spring_len = vc.compute_spring_length(theta)
        spring_length_rel = spring_len - vc.spring_length_natural
        minimum_margin = np.min([np.min(spring_length_rel), minimum_margin])
    return minimum_margin


def constraint_spring_maximum(parm):
    # inequality type constraint. Valid when return > 0
    # compare each spring's maximum reached (delta) length to its maximum rated length
    margin = 1
    for i in range(spring_qty):
        offset = parms_per_spring * i
        vc = VirtualCounterweight(lever_radius=parm[offset+0],
                spring_origin=np.array([parm[offset+1], parm[offset+2], 0]),
                spring_cable_length=parm[offset+3],
                spring_stiffness=spring_stiffness,
                spring_length_natural=spring_length_natural,
                spring_length_max=spring_length_max)
        spring_len = vc.compute_spring_length(theta)
        spring_length_rel = vc.spring_length_max - spring_len
        margin = np.min([np.min(spring_length_rel), margin])
    return margin

constraints = ({"type":"ineq", "fun":constraint_spring_minimum},
                {"type":"ineq", "fun":constraint_spring_maximum})

# parm_guess_unit = [0.025, 0.010, -0.150, 0.020]
# parm_bounds_unit = []

# VC: arm_len, spring_x, spring_y, cable_len
parm_guess = [0.025, -0.030, -0.150, 0.010,   0.025, 0.050, -0.150, 0.010]
# parm_bounds = ([0.015, -0.130, -0.270, 0.000,  0.015, -0.130, -0.270, 0.000], [0.070, 0.130, -0.110, 0.300,   0.070, 0.130, -0.110, 0.300])

parm_bounds = ((0.015, 0.090), (-0.130, 0.130), (-0.270, -0.110), (0.000, 0.300),
               (0.015, 0.090), (-0.130, 0.130), (-0.270, -0.110), (0.000, 0.300))
soln = optim.minimize(residual, x0=parm_guess, bounds=parm_bounds, constraints=constraints)
print(soln.x)


# show the solution

parm = soln.x

counterweights = []
for i in range(spring_qty):
    offset = parms_per_spring * i
    counterweights.append(VirtualCounterweight(lever_radius=parm[offset+0],
        spring_origin=np.array([parm[offset+1], parm[offset+2], 0]),
        spring_cable_length=parm[offset+3],
        spring_stiffness=spring_stiffness,
        spring_length_natural=0.100, # TODO
        spring_length_max=0.200)) # TODO
net_torque = gravity_torque


i = 0
for vc in counterweights:
    tau = vc.compute_torque(theta)
    net_torque += tau
    plt.plot(theta, -tau, label=f"-spring {i}")
    i += 1
plt.plot(theta, net_torque, label="net after compensation")
plt.legend()

net_torque_worst = np.max(np.abs(net_torque))
print(f"  Worst case torque error {net_torque_worst:.2f} Nm is {100*net_torque_worst/gravity_torque_max:.0f}% of the original max")

plt.figure(2)
for i in range(spring_qty):
    offset = parms_per_spring * i
    plt.plot([0, parm[offset + 0], parm[offset + 1]], [0, 0, parm[offset + 2]])
plt.axis("equal")
plt.grid(True)

plt.show()
