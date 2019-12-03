Suit of MATLAB files for analysis, parameter extraction and simulation of
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
Paul trap functions (in particular, extract experimental parameters by fitting to
a series truncation of Mathieu function).


4- simulation
--------------------
Contains code to simulate ODEs and SDEs in MATLAB. Designed to easily
import particle, trap, etc. data. Also contains simulation of LQG (Kalman + 
LQR); in particular for a harmonic oscillator trap.
