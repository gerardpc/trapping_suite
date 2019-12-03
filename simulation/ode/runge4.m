%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Four-order Runge-Kutta method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve an ODE y' = f(t,y), y(t0) = y0
% Used to calculate first steps in Adams predictor-corrector method
% 
% inputs: f_edo -- the function f(t,y), as an inline
%         t_interval -- the interval [t0,tn]
%         y0 -- the initial value
%         h -- the step size
% output: t, y -- the node and the value of y
% 
% GP Conangla, 12/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y] = runge4(f_edo, t_interval, y0, h)

t = t_interval(1):h:t_interval(2);
y = zeros(length(y0), length(t));
y(:,1) = y0;

for n = 1:(length(t) - 1)
    k1 = feval(f_edo, t(n), y(:,n));
    k2 = feval(f_edo, t(n) + h/2, y(:, n) + h/2*k1);
    k3 = feval(f_edo, t(n) + h/2, y(:, n) + h/2*k2);
    k4 = feval(f_edo, t(n) + h, y(:, n) + h*k3);
    y(:, n+1) = y(:, n) + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end
t = t';
y = y';
end
