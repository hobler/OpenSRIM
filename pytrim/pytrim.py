"""PyTRIM aims to be a Python implementation of TRIM.

TRIM (Transport of Ions in Matter) is a widely used software package 
for simulating the interaction of ions with matter, particularly for 
ion implantation in semiconductors. It comes as part of SRIM, see
www.srim.org.

PyTRIM seeks to replicate the core functionalities of TRIM using 
Python, making it more accessible and easier to integrate with other
Python-based tools and workflows.

Currently, the input parameters are hardcoded in this script, but future
versions may include a more user-friendly interface for specifying
simulation parameters. Also, recoils are not yet followed, and only the
mean and the straggling of the penetration depth of the primary ions are
recorded.
"""
from math import sqrt
import time
import numpy as np
import select_recoil
import scatter
import estop
import geometry
import trajectory


start_time = time.time()

nion = 1000             # number of projectiles to simulate

zmin = 0.0              # minimum z coordinate of the target (A)
zmax = 4000.0           # maximum z coordinate of the target (A)
z1 = 5                  # atomic number of projectile
m1 = 11.009             # mass of projectile (amu)
z2 = 14                 # atomic number of target
m2 = 28.086             # mass of target atom (amu)
density = 0.04994       # target density (atoms/A^3)
corr_lindhard = 1.5     # Correction factor to Lindhard stopping power

select_recoil.setup(density)
scatter.setup(z1, m1, z2, m2)
estop.setup(corr_lindhard, z1, m1, z2, density)
geometry.setup(zmin, zmax)
trajectory.setup()

# Initial conditions of the projectile
e_init = 50000.0                         # energy (eV)
pos_init = np.array([0.0, 0.0, 0.0])     # position (A)
dir_init = np.array([0.0, 0.0, 1.0])     # direction (unit vector)

count_inside = 0
mean_z = 0.0
std_z = 0.0
for _ in range(nion):
    pos, dir, e, is_inside = trajectory.trajectory(pos_init, dir_init, e_init)
    if is_inside:
        count_inside += 1
        mean_z += pos[2]
        std_z += pos[2]**2

mean_z /= count_inside
std_z = sqrt(std_z/count_inside - mean_z**2)
print(f"Number of ions stopped inside the target: {count_inside} / {nion}")
print(f"Mean penetration depth of ions stopped inside the target: "
      f"{mean_z:.2f} A")
print(f"Standard deviation of penetration depth: {std_z:.2f} A")

end_time = time.time()
print(f"Simulation time: {end_time - start_time:.2f} seconds")