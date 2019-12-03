%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit_lorentzian(pxx, ww, varargin)
% 
% Fit data to Lorentzian distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%        pxx      -- PSD data
%        ww       -- associated frequency (in omega)
%        varargin -- If 'with_offset', add offset possibility to fit
%
% OUTPUTS: 
%        fit_params     -- Obtained parameters from fit
%        conf_intervals -- Confidence intervals (1 sigma) for parameters
%        ww_fit         -- Omega vector for fit
%        curve_fit      -- Fit vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 7.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[fit_params, conf_intervals, ww_fit, curve_fit] = fit_lorentzian(pxx, ww, varargin)

% %% Load optical tweezer parameters for comparison -> ignore this, it is
% for simulations.
% opt = load_opt_tweezer();

%% Define function to fit: Lorentzian
% x(1) = constant factor (K_factor)
% x(2) = w_0 (resonance frequency)
% x(3) = G_norm factor
% x(4) = noise floor
lorentzian = @(x, w) x(1)./((x(2).^2 - w.^2).^2 + x(3).^2*w.^2) + x(4);

%% Initial guesses for fit parameters
t_ini = toc;
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
x0_ini_guess = [K_factor, w_0, G_norm, 0];

% lower and upper bounds for fit parameters
x0_lb = [K_factor/2, w_0*0.8, G_norm/3, 0];
x0_ub = [K_factor*2, w_0*1.2, G_norm*3, 0];

% Consider offset or not?
offset_or_not = find(ismember(varargin, 'with_offset'));
if(~isempty(offset_or_not))    
    noise_floor = mean(power(freq > 3e6));
    x0_ini_guess = [K_factor, w_0, Q, noise_floor];
    x0_lb = [K_factor/2, w_0*0.8, G_norm/3, 0.8*noise_floor];
    x0_ub = [K_factor*2, w_0*1.2, G_norm*3, 1.2*noise_floor];
end

%% Least squares fit
% some random fitting routine options
options = optimset('Algorithm', 'trust-region-reflective', 'TolFun', 1e-80, 'MaxIter', 2000,'MaxFunEvals',1e4, 'Display', 'none');

% Behold the fit itself
[fit_params, ~, resid, ~, ~, ~, J] = lsqcurvefit(lorentzian, x0_ini_guess, ww, pxx, x0_lb, x0_ub, options);

% get fit parameters
K_factor = fit_params(1);
w_0 = fit_params(2);
G_norm = fit_params(3);

% Get confidence intervals (68%) of obtained parameters
conf_intervals = nlparci(fit_params, resid, 'jacobian', J, 'alpha', 0.32);

%% Print the obtained parameters on command line
expected_K_factor = opt.sigma^2/opt.m^2;
fprintf('Relative error of fit parameters:\n\n');
fprintf('Resonance freq f_0         : %.4e\n', abs(w_0 - opt.w_0)/opt.w_0);
fprintf('G_norm/expected value      : %.4e\n', abs(G_norm - opt.G_norm)/opt.G_norm);
fprintf('K_factor/expected value    : %.4e\n', abs(K_factor - expected_K_factor)/expected_K_factor);

%% Calculate fit (output)
ww_fit = linspace(0, 3*w_0, 5000); % only calculate around resonance
curve_fit =  lorentzian(fit_params, ww_fit);
% curve_fit_2 =  lorentzian([expected_K_factor, opt.freq*2*pi, opt.G_norm, 0], ww_fit);

end












