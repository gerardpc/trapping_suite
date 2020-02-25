%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This serves as a 2nd order function for a numerical method for SDEs
% 1st function: x' = f(1)
% 2nd function: v' = f(2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[f] = diffusion_2n_sde(t, y, varargin)
% use 'understandable' names
x = y(1);
v = y(2);

% parameters opt. tweezer
trap = varargin{1};

% % Definitions: experimental parameters (like in pdf)
% m = 9.2e-18;
% gamma = 3.5e-11;
% w = 2*pi*2e4;
% eps = 6.3e-9;
% T = 295;
% k_B = 1.38065e-23;
% sigma = sqrt(2*k_B*T*gamma);

% No noise in x
f(1) = 0; 
% Stochastic force optical tweezer
f(2) = trap.sigma/trap.m;
% % Stochastic force Paul trap
% f(2) = sigma/m;

f = f';
end