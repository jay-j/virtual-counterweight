#!/usr/bin/env python3
import numpy as np
import scipy.optimize as optim
import matplotlib.pyplot as plt

g = 9.81 # acceleration of gravity, m/s^2
mass = 0.30 # kg, mass of the end effector
L = 0.300   # m, length from pivot to the 

# angle of the main arm, radians
# measured with positive DOWN, straight horizontal zero
theta = np.linspace(-np.radians(70), np.radians(70), 128)
theta_label = "Arm Angle, theta (radians)"

# uncompensated torque to support the arm, Nm
gravity_torque = mass*g * L * np.cos(theta)


# guess starting cam radius
cam_radius3 = 0.040 # m

# calculate the force for the counterweight torque at theta=0
force3 = mass * g * L / cam_radius3
print(f"Force3 = {force3:.1f} N")

# guess a point to afix the other end of the tension spring
spring_origin = np.array([cam_radius3, -0.100, 0])
spring_stiffness = 1000 # N/m
spring_length_natural = 0.05 # m

spring_length3 = force3 / spring_stiffness + spring_length_natural
print(f"Spring_length3 = {spring_length3:.3f} m")

#######################################################################################

def rotate_z(vector, angle_rad):
    vector = vector.reshape((-1,3))
    result = np.zeros_like(vector)
    result[:,0] = vector[:,0]*np.cos(angle_rad) - vector[:,1]*np.sin(angle_rad)
    result[:,1] = vector[:,0]*np.sin(angle_rad) + vector[:,1]*np.cos(angle_rad)
    result[:,2] = vector[:,2]
    return result
    

# cam parametrized by two bezier curves, 8 control points (7 unique)
# A Cubic Bezier Curve; 4 control points; N=3
class Bezier:
    def __init__(self, P0, P1, P2, P3):
        self.N = 3
        self.P0 = P0
        self.P1 = P1
        self.P2 = P2
        self.P3 = P3
    
        # compute control points for the derivative
        self.V0 = self.N * (self.P1 - self.P0)
        self.V1 = self.N * (self.P2 - self.P1)
        self.V2 = self.N * (self.P3 - self.P2)

    def pos(self, t):
        try: # TODO performance! :(
            result = np.zeros((t.shape[0], 3))
        except IndexError:
            result = np.zeros((1, 3))

        for i in range(2):
            result[:,i] = np.power(1-t,3)*self.P0[i] + 3*np.power(1-t,2)*t*self.P1[i] + 3*(1-t)*np.power(t,2)*self.P2[i] + np.power(t,3)*self.P3[i]
        return result

    def slope(self, t):
        # return: [dx/dt, dy/dt]
        # function of Pi, t^2
        try: # TODO performance! :(
            result = np.zeros((t.shape[0], 3))
        except IndexError:
            result = np.zeros((1, 3))

        for i in range(2):
            result[:,i] = np.power(1-t,2)*self.V0[i] + 2*(1-t)*t*self.V1[i] + np.power(t,2)*self.V2[i]
        return result


class Cam:
    def __init__(self, P0, P1, P2, P3, P4, P5, P6):
        self.surface = [Bezier(P0, P1, P2, P3), Bezier(P3, P4, P5, P6)]

    def pos(self, t):
        # trying to be numpy array friendly
        cam1 = self.surface[0].pos(t)
        cam2 = self.surface[1].pos(t-1) # TODO t-1
        combo = np.vstack([cam1[t <= 1],cam2[t > 1]])
        return combo

    def slope(self, t):
        cam1 = self.surface[0].slope(t)
        cam2 = self.surface[1].slope(t-1)
        combo = np.vstack([cam1[t <= 1], cam2[t > 1]])
        return combo


# assume P0 is on the 180. then all you need is a 
P0 = np.array([0, 0.3*cam_radius3, 0]) # TODO guess the radius
P1 = P0 + np.array([0.020, 0.000, 0]) # TODO a fully guessed parameter


P3 = np.array([cam_radius3, -0.030, 0]) # assume this point is on the X axis, radius is variable
P2 = P3 + np.array([0, 0.010, 0]) # TODO a fully guessed parameter
P4 = P3 - np.array([0, 0.010, 0]) # TODO a fully guessed parameter

P6 = np.array([0, -1.1*cam_radius3, 0]) # TODO guess the radius
P5 = P6 + np.array([0.030, 0, 0]) # TODO a fully guessed parameter

cam = Cam(P0, P1, P2, P3, P4, P5, P6)


t = np.linspace(0, 2, 128)
shape = cam.pos(t)
print(shape.shape)


plt.figure(2)
plt.plot(shape[:,0], shape[:,1])
plt.plot(spring_origin[0], spring_origin[1], 'o')


plt.grid(True)
plt.axis("equal")


cam_contact_t = np.zeros(theta.shape) # t value where string touches cam (index match theta)
cam_contact_pos = np.zeros((theta.shape[0], 3)) # cartesian position where string touches cam (index match theta), global position!



def find_tangent_contact(t_guess, cam, spring_origin, th):
    # compute point on cam. compute unit vector from cam point to spring origin
    # compute slope; unit vector representing tangent direction of the curve
    # cross product the two vectors, look for zero using binary search
    # cam = kwargs["cam"]
    # spring_origin = kwargs["spring_origin"]
    # th = kwargs["th"]

    cam_pos = cam.pos(np.array(t_guess))
    cam_pos = rotate_z(cam_pos, th)
    cam_slope = cam.slope(np.array(t_guess))
    cam_slope = rotate_z(cam_slope, th)
    cam_slope /= np.linalg.norm(cam_slope)

    spring_slope = spring_origin - cam_pos
    spring_slope /= np.linalg.norm(spring_slope)

    error = np.cross(spring_slope, cam_slope).reshape(3)[2]
    return error


# TODO this can have a problem where they jump to the other side of the solution space
# solve for the contact conditions, working forwards
index_theta_zero = np.power(theta,2).argmin()
print("theta index..", index_theta_zero)
print("Now starting cam second half")
t_prev = 1.0
for i in range(index_theta_zero, theta.shape[0]):
    th = theta[i]
    # rotate the cam points according to theta
    # search for t where the tangent projection touches the surface
    # warm start based on previous point/guess
    t_prev_pull = 2 - (2 - t_prev)*0.5  # halfway between t_prev and the upper bound, 2
    soln, status = optim.newton(find_tangent_contact, x0=t_prev, x1=t_prev_pull, args=(cam, spring_origin, th), full_output=True, disp=False)
    if not status.converged:
        soln = -1
    cam_contact_t[i] = soln
    cam_contact_pos[i,:] = rotate_z(cam.pos(cam_contact_t[i]), th)

print("Now starting cam first half")
t_prev = 1.0
for i in range(index_theta_zero, -1, -1):
    th = theta[i]
    # rotate the cam points according to theta
    # search for t where the tangent projection touches the surface
    # warm start based on previous point/guess
    soln, status = optim.newton(find_tangent_contact, x0=t_prev, x1=t_prev*0.5, args=(cam, spring_origin, th), full_output=True, disp=False)
    if not status.converged:
        soln = -1
    cam_contact_t[i] = soln
    cam_contact_pos[i,:] = rotate_z(cam.pos(cam_contact_t[i]), th)

print(cam_contact_t)
print(cam_contact_pos.shape)

i = 127
cam_rot = rotate_z(shape, theta[i])
plt.plot(cam_rot[:,0], cam_rot[:,1])

plt.plot([cam_contact_pos[i,0], spring_origin[0]], [cam_contact_pos[i,1], spring_origin[1]])

# plot the first and last contact points
local_start = rotate_z(cam_contact_pos[0,:], -theta[0])
local_end = rotate_z(cam_contact_pos[-1,:], -theta[-1])
plt.plot([local_start[0,0], local_end[0, 0]], [local_start[0,1], local_end[0,1]], 'x')

# calculate update curve length along the entire cam, match index to theta (all that matters is relative length to neutral)
cam_length = np.zeros_like(theta)
for i in range(1, theta.shape[0]):
    # Integrate from the previous step. From the previous stage, add up dY/dt and dX/dt
    t = cam_contact_t[i]
    dt = t - cam_contact_t[i-1]
    cam_length[i] = cam_length[i-1] + np.linalg.norm(cam.slope(t))*dt
    
# offset so cam length=0 at the neutral/zero position
cam_length -= cam_length[index_theta_zero]

free_length = np.zeros_like(theta)
for i in range(theta.shape[0]):
    free_length[i] = np.linalg.norm(cam_contact_pos[i,:] - spring_origin)
free_length -= free_length[index_theta_zero]

plt.figure(3)
plt.plot(theta, cam_length*1e3, label="cam surface")
plt.plot(theta, free_length*1e3, label="free space")
plt.plot(theta, 1e3*(cam_length + free_length), label="net length")
plt.xlabel("theta angle, rad")
plt.ylabel("length, mm")
plt.grid(True)
plt.legend()

spring_length = spring_length3 + cam_length + free_length

# calculate spring force TODO
# F = k*x

# calculate spring torque on arm TODO 
# r cross F



plt.figure(1)
plt.plot(theta, gravity_torque)
plt.xlabel(theta_label)
plt.ylabel("Torque, Nm")
plt.grid(True)

plt.show()
