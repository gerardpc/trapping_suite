%
% rad = particle radius in nm
% std = standard deviation of radius in nm
% material = silica or diamond

function [mass, std_percentage] = mass_particle(rad, std_rad, unit_rad, material, unit_output)

%% Densities and constants
kg_to_amu = 6.022141*10^26; % amus in kg
switch unit_rad
    case 'um'
        rad = rad*1000;
        std_rad = std_rad*1000;
    case 'nm'
        % do nothing
    case 'm'
        rad = rad*10^9;
        std_rad = std_rad*10^9;
end
switch material
    case 'silica'
        density = 2196; % kg/m3
    case 'diamond'
        density = 3500; % kg/m3
    case 'polystyrene'
        density = 1000; % kg/m3
end

%% Rad -> Volume -> mass
v = 4/3*pi*rad^3; % volume in nm3
std_v = 4*pi*rad^2*std_rad;
std_percentage = std_v/v*100;

mass = v/(10^9)^3*density;

%% Print results
fprintf('----------------------------------\n');
switch unit_output
    case 'kg'
        % already calculated
        fprintf('Mass of %s nanoparticle: %.4g +/- %.3g %% kg\n', material, mass, std_percentage);
    case 'amu'
        mass = mass*kg_to_amu;
        fprintf('Mass of %s nanoparticle: %.4g +/- %.3g %% amu\n', material, mass, std_percentage);
end
fprintf('\n');

end

