# Tracing-the-path-of-satellite-in-stable-orbit-by-simulating-in-Python-

Tracing the path of satellite in stable orbit by simulating in
Python
Nipun Aggarwal, Shubham Jadhav
January 2022
1 Introduction
From our childhood almost everyone curious
about our solar system, planets, stars etc. In this
project we did work to explore the same childhood
curiosity of almost every child.In this project we
plotted stable orbits of all planets present in our
solar system by using some machine learning algorithms, We have solved differential equations
for orbits of each planets by using data provided
,and plotted it on a 2-D graph by using Python 3
libraries.
2 Study Areas
Numerical methods to solve differential equations
(in this project finite difference method), Visualizing planets trajectories using python, Some basic classical mechanics concepts
3 Motion under the influence of Central Forces
Central forces are forces between two particles
that depend only on the distance between the
particles and point from one particle to another.
Central forces are conservative . Furthermore,
if we look at a central force between two particles we can separate the motion of the canter
of mass from the relative motion of the particles
reducing the complexity of the problem, in this
case centre(sun) and other planets moving about
their centre of mass, and forming their stable orbits. With mutual forces between planets providing them centripetal acceleration and angular
momentum, and balance of forces due to all planets on a certain palnet results forming an stable
orbit around the sun.
Conservation of angular momentum:
L = r × p (1)
where r is the orbit radius, p is the momentum.
3.1 Effective Potential
The effective potential provides an useful method
for simplifying a three-dimensional central-force
problem down to a one-dimensional.
For a particle of mass m, polar coordinates (r, θ)
L =
1
2
m( ˙r
2 + r
2 ˙θ
2
) − v(r) (2)
is rather interesting.
It involves only the variable r. And it looks a
lot like the equation for a particle moving in one
dimension (labeled by the coordinate r) under the
influence of the potential.
The equations of motion obtained from varying r
and θ are
m ˙r˙ = mr ˙θ
2 − V
′
(r)
d
dt(mr2 ˙θ) = 0
Consider the example where V (r) = Ar2
. This
is the potential for a spring with relaxed length
zero. Then,
Vef f (r) = L
2
2mr2
+ V (r)
Fef f (r) = L
2
mr3
− V
′
(r)
Consider the example where V (r) = Ar2
. This
is the potential for a spring with relaxed length
zero. Then,
Vef f (r) = L
2
2mr2
+ Ar2
(3)
1
Figure 1: Plot of Vef f
The gravitational potential energy of the
earth–sun system is
V (r) = −
α
r
(4)
where α = GMm M = mass of Sun
m = mass of Earth
In the present treatment, we’ll consider the sun
to be bolted down at the origin of our coordinate
system. Since M >> m, this is approximately
true for the Earth–Sun system.
(
1
r
2
dr
dθ )
2 =
2mE
L2
−
1
r
2
+
2mα
rL2
(5)
With all the 1/r terms floating around, it might
be easier to solve for 1/r instead of r. Using
(
dy
dθ )
2 = −
dr
dθ
r
2
(6)
square on the right-hand side gives
(
dy
dθ )
2 = −(y −
mα
L2
)
2 +
2mE
L2
+ (mα
L2
)
2
(7)
If we take z = y −
mα
L2
(
dz
dθ )
2 = −z
2+(mα
L2
)
2
(1+2EL2
mα2
) = −z
2+B
2
(8)
where,
B = (mα
L2
)
r
1 +
2EL2
mα2
(9)
If we take z = Bcos(θ − θ0)
Z
dz
√
B2 − z
2
=
Z
dθ (10)
We know that z =
1
r −
mα
L2
1
r
=
mα
L2
(1 + ϵcosθ) (11)
where
ϵ =
r
1 +
2EL2
mα2
(12)
Equation (11) will be maximum when cosθ = 1
which is when r will be min
rmin =
L
2
mα(1 + ϵ)
(13)
If ϵ < 1
rmax =
L
2
mα(1 − ϵ)
(14)
4 Finite Difference Method
To solve the ODE boundary value problems is the
finite difference method, where we can use finite
difference formulas at evenly spaced grid points
to approximate the differential equations. This
way, we can transform a differential equation into
a system of algebraic equations to solve.
In the finite difference method, the derivatives in
the differential equation are approximated using
the finite difference formulas. We can divide the
the interval of [a,b] into n equal subintervals of
length h as shown in the following figure.
Commonly, we usually use the central difference formulas in the finite difference methods due
to the fact that they yield better accuracy. The
differential equation is enforced only at the grid
points, and the first and second derivatives are:
dy
dx =
yi+1 − yi−1
2h
d
2
dx2
=
yi+1 − 2yi + yi−1
h
2
These finite difference expressions are used to replace the derivatives of y in the differential equation which leads to a system of n+1 linear algebraic equations if the differential equation is linear. If the differential equation is nonlinear, the
algebraic equations will also be nonlinear.
2
Figure 2: Representation of Finite Difference
Method
5 Concepts used for simulation
Newtons law of gravitation and motion : by using
newtons law of gravitation we found out gravitational force on each planet due to every planet
a =
GM
r
2
(15)
By solving above differential equation we found
out coordinates of trajectory of orbit of each
planet and plotted it on a graph by using matplotlib.
6 Code Explanation
Let’s come to the core of the project i.e. the
code. We’ve used the Python programming
language to implement the concepts and theory.
For the sake of simplicity, we first plotted the
orbit for a single planet and then extended the
concept to multiple planets by making a list of
variables instead of a single variable.
As discussed above, we implemented the finite
difference numerical method so as to get the
trajectory of the planets and took the initial parameters from the table for our calculation: Now,
since we have the initial variables x0, y0, vx0,
and vy0, we calculate the desired variables
x1, y1, vx1, vy1 using finite difference method
on the aforementioned formulas generating the
following equations:
x1 = x0 + vx0 × δt
y1 = y0 + vy0 × δt
vx1 = vx0 + F x × δt
vy1 = vy0 + F y × δt
where,
Fx = −
GMx0
(x
2
0 + y
2
0
)
1.5
Fy = −
GMy0
(x
2
0 + y
2
0
)
1.5
denote the interaction force between sun and the
chosen planet.
After calculating the values x1, y1, vx1, and vy1,
we plot these x and y coordinates simultaneously using the matplotlib library. To obtain the
whole trajectory, we assign the obtained variables
x1, y1, vx1, vy1 back to x0, y0, vx0 and vy0 and iteratively repeat the above process until the desired
number of frames. We have successfully obtained
our result for a single planet.
Now, let’s generalise the above code for n planets.
We need to make only minor changes. Firstly,
we’ll make a list of variables x0, y0, vx0, vy0, M, R
instead of maintaining a single variable. Now,
we’ll execute the exact same logic as above, just
for n different planets, one at a time.
So, the above equations would look like:
x1[i
thplanet] = x0[i
thplanet] + vx0[ithplanet] × δt
y1[i
thplanet] = y0[i
thplanet] + vy0[i
thplanet] × δt
vx1[i
thplanet] = vx0[i
thplanet] + F x[i
thplanet] × δt
vy1[i
thplanet] = vy0[i
thplanet] + F y[i
thplanet] × δt
where
F x[i
thplanet] = −Σ
GMδx
(δx2 + δy2
)
1.5
F y[i
thplanet] = −Σ
GMδy
(δx2 + δy2
)
1.5
because we have to consider the interaction between ith planet and all other planets.
Similar to what we have done for a single planet,
3
Figure 3: Initial parameters used for this simulation
Figure 4: Final Simulation results
we plot x1[i
th planet] and y1[i
th planet] simultaneously using the matplotlib library. And we assign the obtained variables x1[i
th planet], y1[i
th
planet], vx1[i
th planet], vy1[i
th planet] back to
x0[i
th planet], y0[i
th planet], vx0[i
th planet] and
vy0[i
th planet] and iteratively repeat the process
for each planet, at different frames until we get
the whole trajectory.
4
import matplotlib.pyplot as plt
# Constants
G = 6.67408 * 1e-11
delta_t = 86400
# Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune
n = 9
x0 = [0, -46000000000, -107480000000, -147095000000, -206620000000, -740520000000, -1352550000000, -2741300000000, -4444450000000]
y0 = [0, 0, 0, 0, 0, 0, 0, 0, 0]
vx0 = [0, 0, 0, 0, 0, 0, 0, 0, 0]
vy0 = [0, -58980, -35260, -30300, -26500, -13720, -10180, -7110, -5500]
M = [1.989*1e30, 0.33011*1e24, 4.8675*1e24, 5.972*1e24, 6.4171*1e23, 1898.19*1e24, 568.34*1e24, 86.813*1e24, 102.413*1e24]
R = [695700000, 2439700, 6051800, 6371000, 3389500, 71492000, 54364000, 24973000, 24341000]
color = ['wo', 'wo', 'wo', 'wo', 'wo', 'mo', 'bo', 'go', 'ro']
frames = 3000
x1 = [0] * n
y1 = [0] * n
vx1 = [0] * n
vy1 = [0] * n
# Solving
for i in range(frames):
for j in range(1, n):
# Calculations
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
# Updations
x0[j] = x1[j]
y0[j] = y1[j]
vx0[j] = vx1[j]
vy0[j] = vy1[j]
# Plotting
plt.plot(x1[j], y1[j], color[j], markersize = 0.5)
# Plotting Sun
plt.plot(0, 0, 'yo', markerSize = 7)
plt.annotate("Sun", (0, 0))
plt.show()
