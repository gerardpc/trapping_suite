%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine for the Labview controler of the Paul trap:
% 1- Read secular motion frequency
% 2- Calculate beta (this and...)
% 3- Get q parameter from beta (...this tell you the stability straightaway)
% 4- "Obtain" (calculate) the number of charges in the particle
%
% Inputs: 
% V := AC voltage amplitude (kV)
% w_driving := AC driving frequency (2*pi*Hz)
% r := Particle radius (nm)
% material := particle material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[q, beta, n] = labview_paultrap(V, w_driving, r, material, varargin)
%% Some preliminary constants and definitions
a = 0; % Mathieu eq. parameter
[m, ~] = mass_particle(r, 0, 'nm', material, 'kg'); % particle mass
%% Calculate secular motion frequency
W_secular = 2*pi*1800;
%% Calculate beta
beta = measure_beta(w_driving, W_secular);
%% Calculate q - Mathieu equation - parameter
q = get_q(a, beta);
%% Number of charges
[~, n] = get_charge(w_driving, m, V, q);
end