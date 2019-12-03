%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_y0()
%
% Pipe (or generate, in case they're random) initial conditions for an SDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%         ini_values -- vector with [x0 v0] OR 'rand_ini', in which
%                           case random y0 is generated
%
% OUTPUTS:
%         y0         -- vector with [x0 v0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 8.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[y0] = get_y0(ini_values)

if(strcmp(ini_values, 'rand_ini'))
    % x0
    x0 = sqrt(1e-11)*normrnd(0, 1);
    % v0
    v0 = sqrt(1.25e-3)*normrnd(0, 1);
    % full vector
    y0 = [x0 v0];
else
    % return input
    y0 = ini_values;
end
end