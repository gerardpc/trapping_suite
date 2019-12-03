%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Omega_z (revolution symmetry axis secular frequency)
%
% w := driving frequency
% m := particle mass
% V := voltage amplitude
% q := Mathieu q_z parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q, n] = get_charge(w, m, V, q)
%% Constants
e = 1.6021766208e-19; % elementary charge
%% characteristic dimension
z0 = 0.7e-3; % this is approx, check real number from simulations
r0 = z0*sqrt(2); % this is the "optimal ratio", check "Charged particle traps" book
d = sqrt(r0^2 + 2*z0^2);
%% Calculations
Q = q*m*d^2*w^2/(4*V); % charge in C
n = Q/e; % number of charges
end

