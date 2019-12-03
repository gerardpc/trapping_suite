%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This serves as a 2nd order function for a numerical method for ODEs
% 1st function: x' = f(1)
% 2nd function: v' = f(2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[f] = test_2nd_order(t, y)
%% use 'understandable' names (DON'T CHANGE THIS SECTION)
x = y(1);
v = y(2);

f(1) = v; 

% SIBA
% eta = 50;
% f(2) = -2*x/(1+ (eta*x^2 + 2)^2);

% Gaussian beam
%% parameters
% general
unit = 1e-9; % nm
h_bar = 1.0546e-34;
c = 3e8;
% laser beam
lambda = 1064*unit; % laser wavelength 
P_0 = 70e-3; % laser power
NA = 0.8; % Gaussian beam NA
% particle 
n = 1.43; % refractive index of particle
a = 136/2*unit; % particle radius
[m_particle, ~] = mass_particle(a, 0, 'm', 'silica', 'kg'); % particle mass
% numerical 
h = 1e-5*a; % discretization size
c_factor = 1.5; % correction factor for high NA in Gaussian beam (set to "1" to eliminate correction)¨

% Gaussian beam and EM parameters definitions
w_0 = lambda/(pi*NA)*c_factor;
z_R = pi*w_0^2/lambda;
eps_0 = 8.854e-12;

%% Field formulas/expressions
% induced bulk dipole
alpha_bulk = 4*pi*a^3*eps_0*(n^2-1)/(n^2+2);

% Gaussian beam
I_0 = P_0*2/(pi*w_0^2);
w = @(z) w_0*sqrt(1 + (z./z_R).^2);
I = @(r, z) I_0*(w_0./w(z)).^2.*exp(-(2*r.^2./w(z).^2));

%grad_I_r = @(r, z) (I(r+h, z) - I(r-h, z));
grad_I_r = @(r, z) (1/12*I(r-2*h, z) - 2/3*I(r-h, z) + 2/3*I(r+h, z) - 1/12*I(r+2*h, z))/h;
grad_I_z = @(r, z) (1/12*I(r, z-2*h) - 2/3*I(r, z-h) + 2/3*I(r, z+h) - 1/12*I(r, z+2*h))/h;

grad_E2_r = @(r, z) grad_I_r(r, z)./(c*eps_0);
grad_E2_z = @(r, z) grad_I_z(r, z)./(c*eps_0);

% dipole force (on a Rayleigh particle (a << lambda) or atom)
F_exact_r = @(r, z, alpha) 1/2*alpha*grad_E2_r(r,z);
F_exact_z = @(r, z, alpha) 1/2*alpha*grad_E2_z(r,z);
f(2) = F_exact_r(x, 0, alpha_bulk)/m_particle;


%% good output format
f = f';
end