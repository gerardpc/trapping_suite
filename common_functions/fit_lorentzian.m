%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit_lorentzian(pxx, ww, varargin)
% 
% Fit data to Lorentzian distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%        pxx      -- PSD data
%        ww       -- associated frequency (in omega)
%
% OUTPUTS: 
%        def_params     -- Obtained parameters from fit
%        conf_intervals -- Confidence intervals (1 sigma) for parameters
%        ww_fit         -- Omega vector for fit
%        curve_fit      -- Fit vector
%        beta           -- Calibration parameter defined as 
%                          2*kB*T*G*w_0^2/K : bits^2/Hz to energy/Hz
%        sigma_beta     -- Uncertainty in beta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 7.2018
% Last version: 23.1.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[def_params, conf_intervals, ww_fit, curve_fit, beta] = fit_lorentzian(ww, pxx)

%% Define function to fit: Lorentzian
% x(1) = constant factor (K_factor)
% x(2) = w_0 (resonance frequency)
% x(3) = G_norm factor
lorentzian = @(x, w) (x(1)./((x(2).^2 - w.^2).^2 + x(3).^2*w.^2));

%% Initial guesses for fit parameters
% smooth pxx to guess initial parameters
s_pxx = smooth(pxx, 1e-3); 
% x(2): guess w_0 as highest point in Sxx
w_0 = ww(s_pxx == max(s_pxx));
% x(1): Estimate K_factor by looking at the mean of low frequencies:
% at low f Sx ~ K_factor/w_0^4
K_factor = mean(s_pxx(ww/(2*pi) > 10e3 & ww/(2*pi) < 30e3))*w_0^4;
% x(3): width of peak
% use the fact that at resonance Sx ~ Q^2/m^2w_0^4
G_norm = sqrt(K_factor/(w_0^2*s_pxx(ww == w_0))); 

%% Least squares fit
% some random fitting routine options
options = optimset('Algorithm', 'trust-region-reflective', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxFunEvals', 1e4, 'Display', 'off');

% Trim data around X mode before the fit
fit_interval = 2*pi*[90e3 150e3];
pxx_trimmed = pxx(ww > fit_interval(1) & ww < fit_interval(2));
ww_trimmed = ww(ww > fit_interval(1) & ww < fit_interval(2));

%% Behold the fit itself
% Fit number 1: K is fixed, optimize w_0 and G_norm
% initial values:
x0_ini_guess = [K_factor, w_0, G_norm];
x0_ini = x0_ini_guess;
x0_lb = [K_factor, w_0*0.8, G_norm/10];
x0_ub = [K_factor, w_0*1.2, G_norm*10];

% Least squares fit 1
[fit_params, ~, ~, ~, ~, ~, ~] = lsqcurvefit(lorentzian, x0_ini, ww_trimmed, pxx_trimmed, x0_lb, x0_ub, options);
%--------------------------------------------------------------------------

% Fit number 2: w_0 is fixed, optimize K and G_norm
% initial values:
x0_ini = [fit_params(1), fit_params(2), fit_params(3)];
x0_lb = [x0_ini(1)/10, fit_params(2), x0_ini(3)/10];
x0_ub = [x0_ini(1)*10, fit_params(2), x0_ini(3)*10];

% Least squares fit 2
[fit_params, ~, ~, ~, ~, ~, ~] = lsqcurvefit(lorentzian, x0_ini, ww_trimmed, pxx_trimmed, x0_lb, x0_ub, options);
%--------------------------------------------------------------------------

% Fit number 3: K is fixed, optimize w_0 and G_norm
% initial values:
x0_ini = [fit_params(1), fit_params(2), fit_params(3)];
x0_lb = [x0_ini(1), fit_params(2)*0.8, x0_ini(3)/10];
x0_ub = [x0_ini(1), fit_params(2)*1.2, x0_ini(3)*10];

% Least squares fit 3
[fit_params, ~, ~, ~, ~, ~, ~] = lsqcurvefit(lorentzian, x0_ini, ww_trimmed, pxx_trimmed, x0_lb, x0_ub, options);
%--------------------------------------------------------------------------

% Fit number 4: w_0 is fixed, optimize K and G_norm
% initial values:
x0_ini = [fit_params(1), fit_params(2), fit_params(3)];
x0_lb = [x0_ini(1)/10, fit_params(2), x0_ini(3)/10];
x0_ub = [x0_ini(1)*10, fit_params(2), x0_ini(3)*10];

% Least squares fit 4
[fit_params, ~, resid, ~, ~, ~, J] = lsqcurvefit(lorentzian, x0_ini, ww_trimmed, pxx_trimmed, x0_lb, x0_ub, options);

% get fit parameters
K_factor = fit_params(1);
w_0 = fit_params(2);
G_norm = fit_params(3);

% Get confidence intervals (68%) for w_0, G_norm
conf_intervals_2 = nlparci(fit_params, resid, 'jacobian', J, 'alpha', 0.32);
conf_intervals_2 = conf_intervals_2(2:end, :); % discard K conf intervals
%--------------------------------------------------------------------------

% Fit number 5: just optimize K, to get conf. intervals
lorentzian_amp = @(x, w) (x(1)./((w_0.^2 - w.^2).^2 + G_norm.^2*w.^2));
% initial values:
x0_ini = K_factor;
x0_lb = K_factor/2;
x0_ub = K_factor*2;

% Least squares fit 5
[fit_params, ~, resid, ~, ~, ~, J] = lsqcurvefit(lorentzian_amp, x0_ini, ww_trimmed, pxx_trimmed, x0_lb, x0_ub, options);

K_factor = fit_params(1);
% Get confidence intervals (68%) for K
conf_intervals_1 = nlparci(fit_params, resid, 'jacobian', J, 'alpha', 0.32);
%--------------------------------------------------------------------------

% Definitive parameter fits
def_params = [K_factor, w_0, G_norm];

% Definitive confidence intervals
conf_intervals = [conf_intervals_1; conf_intervals_2];

%% Calibration
kB = 1.38e-23;
T = 295;
beta = 2*kB*T*G_norm*w_0^2/K_factor; % calibration factor: bits to energy


%% Calculate fit (output)
ww_fit = linspace(0, 3*w_0, 5000); % only calculate around resonance
curve_fit =  lorentzian(def_params, ww_fit);
curve_ini_param =  lorentzian(x0_ini_guess, ww_fit);

figure(101);
clf;
nice_plot(ww, pxx,'','','');
nice_plot(ww_fit, curve_ini_param, '\omega (rad/s)', 'S_x', 'Fit result', 'semilogy', 'color', 'yellow', 'linewidth', '1.5');
nice_plot(ww_fit, curve_fit, '\omega (rad/s)', 'S_x', 'Fit result', 'semilogy', 'color', 'red', 'linewidth', '1.5');

end












