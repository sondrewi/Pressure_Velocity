import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Python script to generate comparison plot between x- and y- velocities found by
# Simple_solve and the results obtained by Ghia et. al. (1982). Given a set of
# files containing the converged solution for Re = 100 and Re = 1000
# this script will generate the comparison plots as seen in the report along the
# vertical and horizontal centre lines

# Assuming
df_Re100 = pd.read_csv('Re100_sols.dat', delim_whitespace=True, header=None)
df_Re100.columns = ['X', 'Y', 'vx', 'vy', 'p']

df_Re1000 = pd.read_csv('Re1000_sols.dat', delim_whitespace=True, header=None)
df_Re1000.columns = ['X', 'Y', 'vx', 'vy', 'p']

# Filter to get rows with Y = 0.05025
u_Re100 = df_Re1[df_Re1['Y'] == 0.05025]
u_Re1000 = df_Re2[df_Re2['Y'] == 0.05025]

# Filter to get rows with X = 0.05025
v_Re100 = df_Re1[df_Re1['X'] == 0.05025]
v_Re1000 = df_Re2[df_Re2['X'] == 0.05025]

# Convert to numpy arrays
u_Re100 = u_Re100.to_numpy()
u_Re1000 = u_Re1000.to_numpy()
v_Re100 = v_Re100.to_numpy()
v_Re1000 = v_Re1000.to_numpy()

# Make arrays of values from Ghia et al, at Re = 100 and Re = 1000. u indicates x-velocity and
# v indicates y-velocity
Ghia_Re100_u = np.array([0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -
                         0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1.0])
Ghia_Re1000_u = np.array([0, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289, -0.27805, -
                          0.10648, -0.06080, 0.05702, 0.18719, 0.33304, 0.46604, 0.51117, 0.57492, 0.65928, 1.0])
Ghia_Re100_v = np.array([0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527,
                         0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0.0])
Ghia_Re1000_v = np.array([0, 0.27485, 0.29012, 0.30353, 0.32627, 0.37095, 0.33075, 0.32235,
                          0.02526, -0.31966, -0.42665, -0.51550, -0.39188, -0.33714, -0.27669, -0.21388, 0.0])
Ghia_y = np.array([0, 0.00547, 0.00625, 0.00703, 0.01016, 0.01719, 0.02813, 0.04531,
                   0.05, 0.06172, 0.07344, 0.08516, 0.09531, 0.09609, 0.09688, 0.09766, .1])
Ghia_x = np.array([0, 0.00625, 0.00703, 0.00781, 0.00938, 0.01563, 0.02266, 0.02344,
                   0.05, 0.08047, 0.08594, 0.09063, 0.09453, 0.09531, 0.09609, 0.09688, .1])

# Plot x-velocities
plt.plot(u_Re100[:, 2], u_Re100[:, 1],
         color='blue', label='Re = 100, Present')
plt.plot(u_Re1000[:, 2], u_Re1000[:, 1],
         color='red', label='Re = 1000, Present')
plt.scatter(Ghia_Re100_u, Ghia_y, marker='^',
            color='blue', label='Re = 100, Ghia et. al.')
plt.scatter(Ghia_Re1000_u, Ghia_y, marker='^', color='red',
            label='Re = 1000, Ghia et. al.')

plt.xlabel('x')
plt.ylabel('y-velocity')
plt.legend()
plt.title('Y-velocity along y=0.05025')
plt.savefig('x_vel_vertical_centre.png')


plt.plot(v_Re100[:, 0], v_Re100[:, 3],
         color='blue', label='Re = 100, Present')
plt.plot(v_Re1000[:, 0], v_Re1000[:, 3],
         color='red', label='Re = 1000, Present')
plt.scatter(Ghia_x, Ghia_Re100_v, marker='^',
            color='blue', label='Re = 100, Ghia et. al.')
plt.scatter(Ghia_x, Ghia_Re1000_v, marker='^', color='red',
            label='Re = 1000, Ghia et. al.')

plt.xlabel('x')
plt.ylabel('y-velocity')
plt.legend()
plt.title('Y-velocity along y=0.05025')
plt.savefig('y_vel_horizontal_centre.png')
