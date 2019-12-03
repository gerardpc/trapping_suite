%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dipole force due to optical gradient in a Gaussian beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No inputs/outputs, generates plots of forces for Gaussian Beam. Recall a
% Gaussian beam assumes the paraxial approximation, which is completely not
% true for large NA (as is the case). There is a correction factor
% in the "parameters" section to approximate better the actual focus.
%
% For reference, check: Teoria, Gerard's notes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last update:
% GP Conangla, 4/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[] = dipole_force()
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
% dipoles
n_dip = 1; % number of dipoles

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

% Gaussian beam Taylor series (quadratic approximation)
I_taylor2 = @(r, z) 2*P_0/pi*(1/w_0).^2 - 8*P_0/pi*(1/w_0)^4*r.^2/2 - 4*P_0*lambda^2/pi^3*(1/w_0)^6*z.^2/2;

grad_I_r_taylor = @(r, z) (1/12*I_taylor2(r-2*h, z) - 2/3*I_taylor2(r-h, z) + 2/3*I_taylor2(r+h, z) - 1/12*I_taylor2(r+2*h, z))/h;
grad_I_z_taylor = @(r, z) (1/12*I_taylor2(r, z-2*h) - 2/3*I_taylor2(r, z-h) + 2/3*I_taylor2(r, z+h) - 1/12*I_taylor2(r, z+2*h))/h;

grad_E2_r_taylor = @(r, z) grad_I_r_taylor(r, z)./(c*eps_0);
grad_E2_z_taylor = @(r, z) grad_I_z_taylor(r, z)./(c*eps_0);

% dipole force (on a Rayleigh particle (a << lambda) or atom)
F_exact = @(r, z, alpha) 1/2*alpha*[grad_E2_r(r,z); grad_E2_z(r,z)];
F_linear = @(r, z, alpha) 1/2*alpha*[grad_E2_r_taylor(r,z); grad_E2_z_taylor(r,z)];

%% Vectors and matrices
r = linspace(-2000*unit, 2000*unit, 1000);
z = linspace(-2000*unit, 2000*unit, 1000);
[R, Z] = meshgrid(r, z);
Int = I(R, Z);

%% Equation of motion
% k_r restoring force
k = -(1/2*F_linear(h, 0, alpha_bulk) - 1/2*F_linear(-h, 0, alpha_bulk))/h;
k_r = k(1);
k = -(1/2*F_linear(0, h, alpha_bulk) - 1/2*F_linear(0, -h, alpha_bulk))/h;
k_z = k(2);
fr_resonance = sqrt(k_r/m_particle)/(2*pi);
fz_resonance = sqrt(k_z/m_particle)/(2*pi);

% print resonance frequencies
fprintf('----------------------------------\n');
fprintf('Resonance frequencies\n');
fprintf('f_r = %.3g  kHz\n', fr_resonance*1e-3);
fprintf('f_z = %.3g  kHz\n', fz_resonance*1e-3);
fprintf('\n');

%% Bulk dipole
% calculate forces in vectors along main axis
% exact
force_r = F_exact(r, 0, alpha_bulk);
force_r = force_r(1,:);
force_z = F_exact(0, z, alpha_bulk);
force_z = force_z(2,:);
% linear
force_r_linear = F_linear(r, 0, alpha_bulk);
force_r_linear = force_r_linear(1,:);
force_z_linear = F_linear(0, z, alpha_bulk);
force_z_linear = force_z_linear(2,:);

%% 2-level system dipole force (non-radiative decay)
% For reference, check "Quantum and Atom optics", chapter 5.8, mainly eq.
% 5.455, and Les Houches ORI notes
% parameters
lambda_laser = 153000*unit; % laser wavelength
lambda_0 = 1532*unit; % atomic resonance
omega_laser = c*2*pi/lambda_laser;
omega_0 = c*2*pi/lambda_0;
detun = omega_laser - omega_0; % detuning
d = 100; % dipole moment of 2 level system
Gamma = d^2*w_0^3/(3*eps_0*pi*h_bar*c^3); % excited state decay rate
gamma_inh_broad = 0; % Set to 0 for homogeneous broadening
gamma_dephasing = Gamma/2 + gamma_inh_broad; % dephasing
E_r = sqrt(I(0,0)/(c*eps_0)); % electric field m+odulus at (r,z) = (0,0). -> this is for single laser field
Omega_rabi = sqrt(2/3)*d*(2*E_r)/h_bar; % Rabi frequency
s = Omega_rabi^2/(Gamma*gamma_dephasing)*1/(1 + detun^2/gamma_dephasing^2); % saturation parameter
alpha_2level = -2*detun*d^2*Gamma/(3*h_bar*Omega_rabi^2*gamma_dephasing)*s/(1 + s);

% dipole force on a 2-level system
% exact
force_r_2level = F_exact(r, 0, alpha_2level);
force_r_2level = force_r_2level(1,:);
force_z_2level = F_exact(0, z, alpha_2level);
force_z_2level = force_z_2level(2,:);

%% Doppler force



%% Plots
% Intensity at focal plane
figure(1);
clf;
box;
hold on;
pcolor(R, Z, Int);
shading interp;
grid off;
xlabel('r (nm)','fontsize', 24, 'FontName', 'CMU Sans Serif');
ylabel('z (nm)','fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
xlim([min(r), max(r)]);
ylim([min(z), max(z)]);

% Intensity along r, z
figure(2);
clf;
% along r
subplot(2,1,1);
box;
hold on;
grid on;
plot(r, I(r, 0),'linewidth',3,'color','b');
plot(r, I_taylor2(r, 0), 'linewidth',3, 'color','r');
xlabel('r (m)','fontsize', 24, 'FontName', 'CMU Sans Serif');
ylabel('Intensity (W/m^2)','fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
ylim([0, inf]);
xlim([min(r), max(r)]);
% along z
subplot(2,1,2);
box;
hold on;
grid on;
plot(z, I(0, z),'linewidth',3, 'color','b');
plot(z, I_taylor2(0, z),'linewidth',3, 'color','r');
xlabel('z (m)','fontsize', 24, 'FontName', 'CMU Sans Serif');
ylabel('Intensity (W/m^2)','fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
ylim([0, inf]);
xlim([min(z), max(z)]);

% Gradient dipole force
figure(3);
clf;
% along r
subplot(2,1,1);
box;
hold on;
grid on;
title('Bulk dipole force');
plot(r, force_r, 'linewidth',3, 'color','b');
plot(r, force_r_linear, 'linewidth',3, 'color','r');
xlabel('r (m)','fontsize', 24, 'FontName', 'CMU Sans Serif');
ylabel('Force (N)','fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
ylim(1.5*[min(force_r), max(force_r)]);
xlim([min(r), max(r)]);
% along z
subplot(2,1,2);
box;
hold on;
grid on;
title('Bulk dipole force');
plot(z, force_z,'linewidth',3, 'color','b');
plot(z, force_z_linear,'linewidth',3, 'color','r');
xlabel('z (m)','fontsize', 24, 'FontName', 'CMU Sans Serif');
ylabel('Force (N)','fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
ylim([0, inf]);
ylim(1.5*[min(force_z), max(force_z)]);
xlim([min(z), max(z)]);

% Gradient dipole force on 2-level system
figure(4);
clf;
% along r
subplot(2,1,1);
box;
hold on;
grid on;
title('2-level dipole force');
plot(r, force_r_2level, 'linewidth',3, 'color','b');
xlabel('r (m)','fontsize', 24, 'FontName', 'CMU Sans Serif');
ylabel('Force (N)','fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
ylim(1.5*[min(force_r_2level), max(force_r_2level)]);
xlim([min(r), max(r)]);
% along z
subplot(2,1,2);
box;
hold on;
grid on;
title('2-level dipole force');
plot(z, force_z_2level,'linewidth',3, 'color','b');
xlabel('z (m)','fontsize', 24, 'FontName', 'CMU Sans Serif');
ylabel('Force (N)','fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
ylim([0, inf]);
ylim(1.5*[min(force_z_2level), max(force_z_2level)]);
xlim([min(z), max(z)]);
end

















