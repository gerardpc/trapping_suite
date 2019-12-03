%
% rad = particle radius in nm
% std = standard deviation of radius in nm
% material = silica or diamond

function [alpha] = polarizability(rad, std_rad, unit_rad, material, unit_output)

%% Densities and constants
alpha = 2*pi*r_part^3/(3*10^8)*(n^2-1)/(n^2+2)*power_term;

end

