%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load_paul_trap()
%
% Read Paul trap parameters from file 'particle_paultrap.txt'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%         paul        -- struct with fields (e.g opt.freq or opt.gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 8.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[paul] = load_paul_trap()
% read file of optical tweezer parameters
file_ID = fopen('particle_paultrap.txt');
text_opt = textscan(file_ID, '%s', 60);
fclose(file_ID);

% put numerical values in vector
values_opt = str2double(text_opt{1});
values_opt(isnan(values_opt) == 1) = [];

% parameters from file...
T = values_opt(1);
k_B = values_opt(2);
r = values_opt(3);
Q = values_opt(4);
pressure = values_opt(5);
density = values_opt(6);
freq = values_opt(7);
V = values_opt(8);
d = values_opt(9);

% ... and calculate rest of parameters
m = 4/3*pi*r^3*density;
gamma_air = 6*pi*r*1.84e-5;
gamma = gamma_air*pressure/1010;
G_norm = gamma/m;
w_0 = 2*pi*freq;
sigma = sqrt(2*k_B*T*gamma);
eps = 1.602176620e-19*Q*V/d^2;
% as defined in izmailov, PRE 1995
beta = 4*Q*V*1.602e-19/(m*(d*w_0)^2);
fprintf('Value of Beta: %.2e\n', beta);

% put values of parameters under struct with understandable names
paul = struct('T', T, 'k_B', k_B, 'r', r,...
    'pressure', pressure, 'density', density, 'freq', freq,...
    'm', m, 'gamma_air', gamma_air, 'gamma', gamma, 'G_norm', G_norm, 'w_0', w_0, ...
    'eps', eps, 'sigma', sigma);
end