%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HARMONIC OSCILLATOR FUNCTION:
% get parameters (check outputs to see exactly what)
% and response function to drivings at different frequencies
%
% Equation of h.o. is defined as:
%
% mx'' + m*Gamma*x' + m*w_0^2*x = sigma_1*eta(t),   (1)
%
% Where eta(t) is a white noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%       m          := mass of particle
%       Q          := quality factor of oscillator
%       w_0        := natural frequency %
% varargin:
%       'radius'   := radius in m (e.g. 'radius', '100e-6')
%       'pressure' := pressure in mbar (e.g. 'pressure', '1000')
%
% OUTPUTS: 
%       harm_parameters    := for cleanliness, all the following parameters
%                             are inside this struct:
%         m                  := mass of harmonic oscillator (maybe has been
%                               given as input, maybe has been calculated 
%                               from spherical particle)
%         Gamma              := damping as defined in eq. (1) above
%         w_r                := resonance frequency (omega of resp_max)
%         resp_max           := value of max. of power transfer function
%         Q                  := quality factor (maybe already given as input)
%         HWHM               := at high Q, HWHM of resonance
%         E_X2               := Expected energy by equipartition (integrated
%                               from driving by Brownian noise)
%         max_force_sens_ASD := ASD of optimum force sensitivity 
%         max_acc_sens_ASD   := ASD of acceleration sensitivity (not taking
%                               into account measurement noise)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMENTS: 
% The quality of the oscillator is defined by the parameter Q.
% However, this parameter can be overwritten by specifying the radius of
% our particle and the pressure. Then, by using Stoke's formula the damping
% Gamma can be estimated. If this is the case, the inputs "m" and "Q" will
% not be used, and instead will be recalculated.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, last update 13.11.2018
%
% Reference: Gerard's theory notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [harm_parameters] = harm_oscillator_response(m, Q, w_0, varargin)
%% preliminary definitions
Gamma = w_0/Q;
kB = 1.3806485e-23; % Boltzmann constant
T = 20 + 273.15; % Temperature of bath

%% variable arguments in
% pressure accelerometer
pressure_pos = find(ismember(varargin, 'pressure'));
if(~isempty(pressure_pos))
    pressure = str2double(varargin{pressure_pos + 1});
end
% radius particle for accelerometer
radius_pos = find(ismember(varargin, 'radius'));
if(~isempty(radius_pos))
    part_radius = str2double(varargin{radius_pos + 1});
    % calculate damping
    density = 1850; % kg/m^3, Stöber silica
    m = 4/3*pi*part_radius^3*density;
    eta = 18.6e-6*pressure/1000; % air viscosity
    gamma_stokes = 6*pi*eta*part_radius; % damping from stokes formula
    Gamma = gamma_stokes/m;
    Q = w_0/Gamma;
end

%% Calculations 
% calculate freq response
f_lim = 1000*w_0;
ww = linspace(0, f_lim, 300000); % omega vector
dw = ww(2) - ww(1);
H2_ptf = @(w, w_0, Q) 1/m^2./((w_0^2 - w.^2).^2 + w.^2*w_0^2/Q^2); % power transfer function
h2 = H2_ptf(ww, w_0, Q); % response

% max response and resonance
resp_max = 4*Q^4/(m^2*(-1 + 4*Q^2)*w_0^4);
w_r = w_0*sqrt(1 - 1/(2*Q^2)); % resonance freq

% bandwidth
HWHM = Gamma/2; % only if Q is high enought (about say > 10)

% noise and sensitivity
sigma_1 = sqrt(2*m*Gamma*kB*T); % thermal noise
sigma_2 = sigma_1*1e30; 
g_min = sqrt(sigma_1^2 + sigma_2^2./H2_ptf(ww, w_0, Q)); % force sensitivity (as root mean square)
max_force_sens_ASD = max(g_min);
max_acc_sens_ASD = sqrt(2*kB*T*Gamma/m);

% integral of Brownian noise response (equipartition check) 
E_X2 = sigma_1^2*trapz(h2)*dw;

% Put everything in a struct that is easy to read for output
harm_parameters = struct('m', m, 'Gamma', Gamma, 'sigma', sigma_1, 'w_r', w_r,...
    'resp_max', resp_max, 'Q', Q, 'HWHM', HWHM, 'E_X2', E_X2,...
    'max_force_sens_ASD', max_force_sens_ASD, 'max_acc_sens_ASD', max_acc_sens_ASD);

%% Plot results
% power transfer function with linear axis
figure(1);
clf;
hold on;
semilogx(ww/(2*pi), h2, 'Color', [1 0 0], 'LineWidth', 3);
ylabel('H(\omega) (m^2/N^2)','fontsize', 24, 'FontName', 'CMU Sans Serif');
title({'Power transfer function'});
xlabel('f (Hz)', 'fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
grid on;
box;
xlim([0 3*w_0/(2*pi)]);

% power transfer function with log axis (bode plot)
figure(2);
clf;
box;
loglog(ww/(2*pi), (h2), 'Color', [1 0.5 0], 'LineWidth', 3);
hold on;
ylabel('H(\omega) (m^2/N^2)','fontsize', 24, 'FontName', 'CMU Sans Serif');
title({'Power transfer function'});
xlabel('f (Hz)', 'fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
grid on;
xlim([0 10*w_0/(2*pi)]);

% response to Brownian noise with log y-axis
figure(3);
% clf;
hold on;
box;
semilogy(ww/(2*pi), sigma_1^2*(h2), 'Color', [0 0.5 0.7], 'LineWidth', 3);
ylabel('$S_x$ (W/Hz)','fontsize', 24, 'FontName', 'CMU Sans Serif', 'Interpreter', 'Latex');
title({'Response to Brownian noise'});
xlabel('f (Hz)', 'fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
grid on;
xlim([0 2*w_0/(2*pi)]);
ylim([10^-24 10^-18]);

% force sensitivity
figure(4);
clf;
hold on;
plot(ww/(2*pi), g_min, 'Color', [0 0 1], 'LineWidth', 3);
title('Force sensitivity (rms)');
ylabel('N/Hz^{1/2}','fontsize', 24, 'FontName', 'CMU Sans Serif');
xlabel('f (Hz)','fontsize', 24, 'FontName', 'CMU Sans Serif');
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif'); % Big labels in axis
grid on;
box;
xlim([0 3*w_0/(2*pi)]);
end






