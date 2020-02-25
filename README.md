Suite of MATLAB files for analysis, parameter extraction and simulation of
data of trapped particles.


1- common_functions
--------------------
Contains functions for (beautiful) plotting, spectral density estimation,
data importing and other commonly used functions.


2- opt_trap
--------------------
Optical trap functions (polarizability, mass, dipole force in Gaussian Beam).


3- paul_trap
--------------------
Paul trap functions (in particular, obtain charge by fitting to
a series truncation of Mathieu function). Main function is called labview_paultrap.m,
which sequentally calls the others.

4- simulation
--------------------
Contains code to simulate ODEs and SDEs in MATLAB. Designed to easily
import particle, trap, etc. data. 
It is divided in three folders:

1. ODE: numerical simulation of noise-less ODE. Main function is full_numerical_ode,
which can be called as:

        >> full_numerical_ode([0, 1e-4], 1e-7, '2nd_order', [0 0]);

The function in the ODE can be changed in the file test_2nd_order.m

2. SDE: numerical simulation of SDE (i.e., ODE + noise). Main function is full_numerical_sde,
which can be called as:

        >> full_numerical_sde([0, 1e-4], 1e-7, '2nd_order', [0 0], 'runge_kutta', 100, 'optical_tweezer','plot')

The function in the ODE can be changed in the files drift_2n_sde.m and diffusion_2n_sde.m

3. Kalman related: Code to simulate Kalman filters, LQR, LQG. 
In particular for a harmonic oscillator trap.

