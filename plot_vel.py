import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.integrate import solve_ivp
from scipy.interpolate import griddata

# Read Solutions file
sols = 'Sols_Re100.dat'
data = np.loadtxt(sols)

# Get values from file
x_coords = data[:, 0]
y_coords = data[:, 1]
x_vel = data[:, 2]
y_vel = data[:, 3]
points = np.vstack((x_coords, y_coords)).T

delta_x = 0.1 / 200
delta_y = 0.1 / 200


# Convert to grid
grid_size = 200
x_grid = x_coords.reshape((200, 200))
y_grid = y_coords.reshape((200, 200))
u_grid = x_vel.reshape((200, 200))
v_grid = y_vel.reshape((200, 200))

# Select every 10th cell, including edges
indices = np.arange(0, grid_size, 8)
indices = np.unique(np.append(indices, [0, grid_size - 1]))

start_points = [[0.07, 0.07]]

# Calculate the norm of velocity
velocity_norm = np.sqrt(u_grid**2 + v_grid**2)

# Function to compute the velocity at any point by fitting a plane from adjacent points


def vel_field(t, xy):
    if (xy[0] < 0 or xy[0] > 0.1) or (xy[1] < 0 or xy[1] > 0.1):
        return [0, 0]

    else:
        cell_x = int(xy[0] / delta_x)
        cell_y = int(xy[1] / delta_y)

        point = np.array([xy[0], xy[1]])
        cell_cent = np.array(
            [(cell_x + 0.5) * delta_x, (cell_y + 0.5) * delta_y])

        up_down = 1 if (cell_cent[1] < point[1]) else -1
        right_left = 1 if (cell_cent[0] < point[0]) else -1

        boundary_x = ((cell_x + right_left == -1)
                      or (cell_x + right_left == grid_size))
        boundary_y = ((cell_y + up_down == -1)
                      or (cell_y + up_down == grid_size))

        boundary_x_coord = 0 if (right_left < 0) else 0.1
        boundary_y_coord = 0 if (up_down < 0) else 0.1

        y_b_value = 0 if (up_down < 0) else 1

        A = np.zeros((4, 3))
        b_u = np.zeros(4)
        b_v = np.zeros(4)

        A[:, 2] = 1

        A[0, 0] = cell_cent[0]
        A[0, 1] = cell_cent[1]

        b_u[0] = u_grid[cell_y, cell_x]
        b_v[0] = v_grid[cell_y, cell_x]

        A[1, 0] = (cell_cent[0] + right_left *
                   delta_x) if (not boundary_x) else (boundary_x_coord)
        A[1, 1] = (cell_cent[1]) if (not boundary_x) else (point[1])
        b_u[1] = u_grid[cell_y, cell_x + right_left] if (not boundary_x) else 0
        b_v[1] = v_grid[cell_y, cell_x + right_left] if (not boundary_x) else 0

        A[2, 0] = (cell_cent[0]) if (not boundary_y) else (point[0])
        A[2, 1] = (cell_cent[1] + up_down *
                   delta_x) if (not boundary_y) else (boundary_y_coord)
        b_u[2] = u_grid[cell_y +
                        up_down, cell_x] if (not boundary_y) else y_b_value
        b_v[2] = v_grid[cell_y + up_down, cell_x] if (not boundary_y) else 0

        if(boundary_x or boundary_y):
            A_inv = np.linalg.inv(A[:3, :])

            u_coefs = A_inv @ b_u[:3]
            v_coefs = A_inv @ b_v[:3]

            u = u_coefs[0] * point[0] + u_coefs[1] * point[1] + u_coefs[2]
            v = v_coefs[0] * point[0] + v_coefs[1] * point[1] + v_coefs[2]

            return [u, v]

        else:
            A[3, 0] = (cell_cent[0] + right_left * delta_x)
            A[3, 1] = (cell_cent[1] + up_down * delta_y)
            b_u[3] = u_grid[cell_y + up_down, cell_x + right_left]
            b_v[3] = v_grid[cell_y + up_down, cell_x + right_left]

            A_pseudo_inv = np.linalg.pinv(A)

            u_coefs = A_pseudo_inv @ b_u
            v_coefs = A_pseudo_inv @ b_v

            u = u_coefs[0] * point[0] + u_coefs[1] * point[1] + u_coefs[2]
            v = v_coefs[0] * point[0] + v_coefs[1] * point[1] + v_coefs[2]

            return [u, v]


plt.figure(figsize=(10, 8))

# Generate streamlines
for start in start_points:
    sol = solve_ivp(vel_field, [0, 25], start,
                    method='RK45', max_step=0.01, vectorized=False)
    # Plot streamline
    plt.plot(sol.y[0], sol.y[1], 'r-', linewidth=0.7)
    print("point done")

# Color plot of the velocity magnitude
norm = Normalize(vmin=velocity_norm.min(), vmax=velocity_norm.max())
plt.pcolormesh(x_grid, y_grid, velocity_norm, shading='auto')
cbar = plt.colorbar(orientation='vertical')
cbar.set_label(r'$|\mathbf{u}|$', rotation=0, labelpad=5, fontsize='large')

u_grid /= velocity_norm
v_grid /= velocity_norm

# quiver plot
plt.quiver(x_grid[indices[:, None], indices], y_grid[indices[:, None], indices],
           u_grid[indices[:, None], indices], v_grid[indices[:, None], indices],
           color='w', scale_units='xy', scale=250, angles='xy')

plt.title('Velocity field: Re = 100')
plt.xlim(0, 0.1)
plt.ylim(0, 0.1)
plt.savefig('vel_Re100.png')
