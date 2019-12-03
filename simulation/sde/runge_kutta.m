%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge - Kutta method for SDE (strong order 1) without derivatives of b
%
% Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%         a          -- the function f(t,y) (string name, as 'sin')
%         b          -- matrix with sigmas of noise term
%         t_interval -- the interval [t0,tn]
%         y0         -- the initial values
%         h          -- the step size
%         varargin   -- If included, contains at varargin{1} the struct
%                    -- opt, from load_opt_tweezer
% OUTPUTS:
%         t, y       -- the node and the value of y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:
% https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_method_%28SDE%29#Variation_of_the_Improved_Euler_is_flexible
%
% GP Conangla, 7/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y] = runge_kutta(a, b, t_interval, y0, h, varargin)

% set t, y vectors and initial values
t = t_interval(1):h:t_interval(2);
y = zeros(length(y0), length(t));
y(:,1) = y0;
dt = h; % nothing special here, just nomenclature
sk_vec = [-1 1];

% Runge - Kutta method
for n = 1:(length(t)-1)
    dW_t = normrnd(0,1,[length(y0),1])*sqrt(dt);
    S_k = sk_vec(randi(2));
    k_1 = feval(a, t(n), y(:,n), varargin{1})*dt + (dW_t - S_k*sqrt(dt)).*feval(b, t(n), y(:,n), varargin{1});
    k_2 = feval(a, t(n + 1), y(:,n) + k_1, varargin{1})*dt + (dW_t + S_k*sqrt(dt)).*feval(b, t(n), y(:,n), varargin{1});
    y(:, n+1) = y(:,n) + 1/2*(k_1 + k_2);
end

% return column vectors with results
t = t';
y = y';
end
