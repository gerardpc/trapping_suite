%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler - Maruyama method for SDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve a langevin type SDE dY = a(t,Y)dt + b*dW_t, Y(t0) = Y0
% inputs: a          -- the function f(t,y) (string name, as 'sin')
%         b          -- matrix with sigmas of noise term
%         t_interval -- the interval [t0,tn]
%         y0         -- the initial values
%         h          -- the step size
% output: t, y       -- the node and the value of y
%
% GP Conangla, 12/2015
%
% Obs (7/2018): superseded by runge_kutta (order 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y] = euler_maruyama(a, b, t_interval, y0, h, varargin)

% set t, y vectors and initial values
t = t_interval(1):h:t_interval(2);
y = zeros(length(y0), length(t));
y(:,1) = y0;
dt = h; % nothing special here, just nomenclature

% Euler - Maruyama method
for n = 1:(length(t)-1)
    dW_t = normrnd(0,1,[length(y0),1])*sqrt(dt);
    y(:, n+1) = y(:,n) + feval(a, t(n), y(:,n), varargin{1})*dt + feval(b, t(n), y(:,n), varargin{1}).*dW_t;
end

% return column vectors with results
t = t';
y = y';
end
