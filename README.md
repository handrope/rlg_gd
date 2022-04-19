# rlg_gd
Simulation of the gradient descent for the Random Lorentz Gas in high dimensions

Code used in the paper https://arxiv.org/abs/2201.01161

Simulation of the gradient descent dynamics of the tracer, Eqs. (45,46)

Compilation:

	g++ -o RLG RLG.cpp -lm
	
Input parameters: d, phi, A

(Dimensions, rescaled density, cutoff radius)

Simulates nrun trajectories with different realizations of obstacles positions.
Trajectories stop when the velocity is lower than vstop and time is higher than tfin

MSD, contacts, pressure and energy vs. time are averaged and printed in "ave.dat"
