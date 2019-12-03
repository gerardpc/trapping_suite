%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This serves as a 1st order function for a numerical method for SDEs
% 1st function: x' = f(1)
% 2nd function: v' = f(2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[f] = drift_1n_sde(t, y, varargin)
% use 'understandable' names
x = y(1);

% parameters opt. tweezer
opt = varargin{1};

% function
f = -x; %any function of t and y;
end