import numpy as np

from solveMotion import solveMotion

rS = np.float64(0.795) # Radius of the sphere (in mm)
R_f = 10 # 52.4/rS # Number of raddi in half a length of the membrane (dimensionless)
Tm = np.float64(27) # np.float64(73) # Tesion of the membrane (mg/ms^2)
g = 9.80665e-3 # Gravity of earth (mm/ms^2)
mu = 1.68e-2 # Density of membrane per unit of area (mg/mm^2) (Sara Wrap membrane)
mS = 73.8 * 4 * np.pi * (rS**3) / 3 # Mass of the ball (mg) (7.8 is the ball's 
# density in mg/mm^3)

solveMotion(rS = rS, Tm = Tm, v_k = -0.50323, \
    plotter = True, recordTime=0.04)

