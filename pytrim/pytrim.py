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
import time
import numpy as np
import select_recoil
import scatter
import estop
import geometry
import cascade
import statistics
import mytypes


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

# Setup modules
select_recoil.setup(density)
scatter.setup(z1, m1, z2, m2)
estop.setup(corr_lindhard, z1, m1, z2, density)
geometry.setup(zmin, zmax)
cascade.setup()
statistics.setup(nspec=1, nbin=40, limits=(0.0, 4000.0))

# Initial conditions of the projectile
proj_init = mytypes.Projectile(
    e = 50000.0,                         # energy (eV)
    pos = np.array([0.0, 0.0, 0.0]),     # position (A)
    dir = np.array([0.0, 0.0, 1.0]),     # direction (unit vector)
)

# Simulate the trajectories
for _ in range(nion):
    proj = proj_init.copy()
    proj_lst, is_inside = cascade.trajectory(proj)
    if is_inside:
        statistics.score(proj_lst[-1])

end_time = time.time()
print(f"Simulation time: {end_time - start_time:.2f} seconds")

# Output the results
statistics.print_results()
statistics.plot_results()
