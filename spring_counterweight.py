#!/usr/bin/env python3
import numpy as np
import scipy.optimize as optim
import matplotlib.pyplot as plt

## DESIGN PARAMETERS #########################################################

g = 9.81 # acceleration of gravity, m/s^2
mass_ee = 0.291 + 0.5 * 0.0724 # kg, mass of the end effector + half cone
L_ee = 0.409   # m, length from pivot to the 
mass_beams = 0.467 # kg

# angle of the main arm, radians
# measured with positive DOWN, straight horizontal zero
theta = np.linspace(-np.radians(53.6), np.radians(43.4), 32)
theta_label = "Arm Angle, theta (radians)"

def gravity_torque_(theta):
    return g*(mass_ee*L_ee + mass_beams*L_ee*0.5) * np.cos(theta)

# uncompensated torque to support the arm, Nm
gravity_torque = gravity_torque_(theta)
gravity_torque_max = np.max(gravity_torque)

##  ##########################################################################

def rotate_z(vector, angle_rad):
    # one vector, many angles
    vector = vector.reshape((-1,3))
    result = np.zeros((angle_rad.shape[0], 3))
    result[:,0] = vector[:,0]*np.cos(angle_rad) - vector[:,1]*np.sin(angle_rad)
    result[:,1] = vector[:,0]*np.sin(angle_rad) + vector[:,1]*np.cos(angle_rad)
    result[:,2] = vector[:,2]
    return result


class Spring:
    def __init__(self, stiffness, length_natural, length_max, od, pn, force_init):
        self.stiffness = stiffness
        self.length_natural = length_natural
        self.length_max = length_max
        self.od = od
        self.pn = pn
        self.force_init = force_init
        

# use a model where the lever parametrization is along +X only. Can rotate the entire system later if needed
class VirtualCounterweight:
    def __init__(self, lever_radius, spring_origin, spring_cable_length, spring):
        self.lever = np.array([lever_radius, 0, 0])
        self.spring_origin = spring_origin
        self.static_length = spring_cable_length
        self.spring = spring
    
    def compute_spring_length(self, theta):
        lever = rotate_z(self.lever, theta)
        spring_vec = self.spring_origin - lever
        spring_length = np.linalg.norm(spring_vec, axis=1).reshape((-1,1))
        return spring_length
    
    def compute_force(self, theta):
        # compute new endpoint given radius
        lever = rotate_z(self.lever, theta)

        spring_vec = self.spring_origin - lever
        
        # compute spring direction
        spring_dir = spring_vec / np.linalg.norm(spring_vec, axis=1).reshape((-1,1))

        # compute spring length
        spring_length_delta = np.linalg.norm(spring_vec, axis=1).reshape((-1,1)) - self.static_length - self.spring.length_natural
        # print("spring length delta", spring_length_delta)

        # compute spring force
        spring_force = (self.spring.force_init + spring_length_delta * self.spring.stiffness) * spring_dir
        # print("spring force:", spring_force)
        return spring_force
            
    def compute_force_magnitude(self, theta):
        force = self.compute_force(theta)
        return np.linalg.norm(force, axis=1)
    
    def compute_torque(self, theta):
        # compute new endpoint given radius
        lever = rotate_z(self.lever, theta)

        # compute spring force
        spring_force = self.compute_force(theta)
        # print("spring force:", spring_force)

        # compute torque = r cross F
        torque = lever[:,0]*spring_force[:,1] + lever[:,1]*spring_force[:,0]
        return torque


# should unique spring stiffnesses per spring be allowed? 
# spring_stiffness = 450
# spring_length_natural = 0.100
# spring_length_max = 0.200
parms_per_spring = 4

def vc_list_from_parameters(parm, spring, spring_qty):
    vc_list = []
    for i in range(spring_qty):
        offset = parms_per_spring * i
        vc_list.append(VirtualCounterweight(lever_radius=parm[offset+0],
            spring_origin=np.array([parm[offset+1], parm[offset+2], 0]),
            spring_cable_length=parm[offset+3],
            spring=spring))
    return vc_list


def compute_net_torque(vc_list):
    net_torque = gravity_torque_(theta)
    for vc in vc_list:
        net_torque += vc.compute_torque(theta)
    return net_torque



def constraint_spring_minimum(parm, spring, spring_qty):
    # inequality type constraint. Valid when return > 0
    # spring delta length must be > 0 (always in tension)
    minimum_margin = 1
    for i in range(spring_qty):
        offset = parms_per_spring * i
        vc = VirtualCounterweight(lever_radius=parm[offset+0],
                spring_origin=np.array([parm[offset+1], parm[offset+2], 0]),
                spring_cable_length=parm[offset+3],
                spring=spring)
        spring_len = vc.compute_spring_length(theta)
        spring_length_rel = spring_len - vc.spring.length_natural
        minimum_margin = np.min([np.min(spring_length_rel), minimum_margin])
    return minimum_margin


def constraint_spring_maximum(parm, spring, spring_qty):
    # inequality type constraint. Valid when return > 0
    # compare each spring's maximum reached (delta) length to its maximum rated length
    margin = 1
    for i in range(spring_qty):
        offset = parms_per_spring * i
        vc = VirtualCounterweight(lever_radius=parm[offset+0],
                spring_origin=np.array([parm[offset+1], parm[offset+2], 0]),
                spring_cable_length=parm[offset+3],
                spring=spring)
        spring_len = vc.compute_spring_length(theta)
        spring_length_rel = vc.spring.length_max - spring_len
        margin = np.min([np.min(spring_length_rel), margin])
    return margin

def residual(parm, spring, spring_qty):
    # the function to optimize!
    # six parameters, dimensional locations of the spring
    # given a constant spring
    vc_list = vc_list_from_parameters(parm, spring, spring_qty)
    tau = compute_net_torque(vc_list)

    # use least squares error
    error = np.sum(np.power(tau, 2))
    return error
    

def spring_optimize(spring, spring_qty):
    # Given a spring, find the optimal arrangement for it
    # returns a solution object from scipy.optimize.minimize
    # parameters: arm_len, spring_x, spring_y, cable_len
    parm_guess = [0.025, -0.030, -0.150, 0.010,   0.025, 0.050, -0.150, 0.010]

    parm_bounds = ((0.015, 0.090), (-0.130, 0.130), (-0.270, -0.110), (0.000, 0.300),
                   (0.015, 0.090), (-0.130, 0.130), (-0.270, -0.110), (0.000, 0.300))

    args = (spring, spring_qty)

    constraints = ({"type":"ineq", "fun":constraint_spring_minimum, "args":args},
                    {"type":"ineq", "fun":constraint_spring_maximum, "args":args})

    soln = optim.minimize(residual, x0=parm_guess, bounds=parm_bounds, constraints=constraints, args=args)
    return soln


def show_solution(spring, spring_qty, solution_parameters):
    # For a given solution, plot torque curves and a schematic of it
    parm = solution_parameters

    vc_list = vc_list_from_parameters(parm, spring, spring_qty)
    net_torque = compute_net_torque(vc_list)

    # plot the torque curves
    plt.figure(1)
    plt.plot(theta, gravity_torque, "k", label="Uncompensated gravity")
    plt.xlabel(theta_label)
    plt.ylabel("Torque, Nm")
    plt.grid(True)

    colors = ["b", "g"]
    for i, vc in enumerate(vc_list):
        tau = vc.compute_torque(theta)
        plt.plot(theta, -tau, colors[i], label=f"-spring {i}")
    plt.plot(theta, net_torque, "r", label="net after compensation")
    plt.legend()

    # performance stat
    net_torque_worst = np.max(np.abs(net_torque))
    print(f"  Worst case torque error {net_torque_worst:.2f} Nm is {100*net_torque_worst/gravity_torque_max:.0f}% of the original max")

    # draw the schematic at theta=0
    plt.figure(2)
    for i in range(spring_qty):
        offset = parms_per_spring * i
        plt.plot([0, parm[offset + 0], parm[offset + 1]], [0, 0, parm[offset + 2]], colors[i])
    plt.axis("equal")
    plt.grid(True)


    # compute force
    plt.figure(3)
    for i, vc in enumerate(vc_list):
        force = vc.compute_force_magnitude(theta)
        plt.plot(theta, force, colors[i], label=f"spring {i}")
    plt.xlabel(theta_label)
    plt.ylabel("Force, N")
    plt.grid(True)
    plt.legend() 

    plt.show()


if __name__ == "__main__":
    spring_qty = 2
    test_spring = Spring(stiffness=450, length_natural=0.100, length_max=0.200, od=0.02, pn="5k3", force_init=20)
    soln = spring_optimize(test_spring, spring_qty)

    show_solution(test_spring, spring_qty, soln.x)

