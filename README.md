Tracing the path of satellite in stable orbit by simulating in Python


#Introduction
This project simulates the stable orbits of planets in our solar system using Python. The project explores central forces, numerical methods, and classical mechanics concepts to trace the paths of planets around the sun. The core objective is to solve differential equations governing planetary motion and visualize the orbits using Python libraries.

#Study Areas

Numerical methods to solve differential equations (finite difference method).
Visualization of planetary trajectories using Python.
Classical mechanics concepts related to central forces and angular momentum.

#Concepts Used : 
Central Forces and Conservation of Angular Momentum
Central forces are conservative and depend only on the distance between particles.
The motion of planets is modeled by separating the center of mass and relative motion.
The effective potential simplifies the three-dimensional central-force problem to one dimension.
Finite Difference Method
Approximates differential equations using finite difference formulas at grid points.
Converts differential equations into algebraic equations for easier solution.
Gravitational Force Calculation
Newton's law of gravitation is used to calculate the force on each planet due to others.
Differential equations are solved to determine the coordinates of planetary orbits.

#Code Explanation :
The code is implemented in Python, utilizing the finite difference method to calculate and plot the orbits of planets. Initially, the orbit of a single planet is plotted, and later the concept is extended to multiple planets.

#Steps Involved:
Initialize Parameters: Initial positions, velocities, and masses of the planets are set based on known data.
Finite Difference Method:
Calculate new positions and velocities using finite difference approximations.
Consider the interaction forces between planets.
Plot Trajectories:
The x and y coordinates are plotted using Matplotlib.
Iteratively repeat the process for all planets and frames.
#Code Snippet

import matplotlib.pyplot as plt

# Constants
G = 6.67408 * 1e-11
delta_t = 86400
n = 9

# Initial positions, velocities, and masses of planets
x0 = [0, -46000000000, -107480000000, ... ]
y0 = [0, 0, 0, ... ]
vx0 = [0, 0, 0, ... ]
vy0 = [0, -58980, -35260, ... ]
M = [1.989*1e30, 0.33011*1e24, ... ]

# Plotting orbits
for i in range(frames):
    for j in range(1, n):
        x1[j] = x0[j] + vx0[j] * delta_t
        y1[j] = y0[j] + vy0[j] * delta_t
        Fx = 0
        Fy = 0
        for k in range(n):
            if(j != k):
                delta_x = x0[j] - x0[k]
                delta_y = y0[j] - y0[k]
                Fx -= G*M[k]*(delta_x) / ((delta_x**2 + delta_y**2) ** 1.5)
                Fy -= G*M[k]*(delta_y) / ((delta_x**2 + delta_y**2) ** 1.5)
        vx1[j] = vx0[j] + Fx * delta_t
        vy1[j] = vy0[j] + Fy * delta_t
        x0[j] = x1[j]
        y0[j] = y1[j]
        vx0[j] = vx1[j]
        vy0[j] = vy1[j]
        plt.plot(x1[j], y1[j], 'wo', markersize=0.5)

# Plotting the Sun
plt.plot(0, 0, 'yo', markersize=7)
plt.annotate("Sun", (0, 0))
plt.show()
Results
The simulation successfully plots the trajectories of all planets in the solar system, showcasing the stable orbits formed due to the balance of gravitational forces.

Conclusion
This project demonstrates how classical mechanics and numerical methods can be combined to simulate and visualize planetary orbits. The finite difference method proves to be an effective tool for solving differential equations governing such complex systems.

Authors
Nipun Aggarwal
Shubham Jadhav
