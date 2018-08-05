# elastic_rod_fracture_simulation
Code accompanying "Controlling fracture cascades through twisting and quenching"

Code 'elastic_rod_fracture.m' runs in MATLAB_R2018a
Simulation code takes the following input parameters (in cgs units) from the file "input_data.txt":

radius, density, Young's modulus, Poisson's ratio, twist damping parameter, minimum fragment size, end-to-end distance at fracture, total twist in rod (degrees)

Simulation proceeds as follows:

  1) End-to-end bending distance is used to choose initial condition from "initial_rod_shapes.mat"
  2) Rod is evolved quasistatically until twist in rod matches input twist
  3) Simulation evolves rod assuming that the maximum initial stress in the rod is the critical stress
  4) Simulation creates video "rod_fracture.avi"
  5) Output array X, of size 3 x 51 x 2001 is stored in "rod_fracture.mat". The first two coordinates of the array are the 51      points of the rod. The last coordinate is the time coordinate.
  
 Plotting tool is a modified version of the "tubeplot" function from the matlab file exchange (See disclaimer in source code, 'elastic_rod_fracture.m')
