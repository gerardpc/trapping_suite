%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Four-step Adams predictor-corrector method 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve an ODE y' = f(t,y), y(t0) = y0
%
% inputs: f_edo      -- the function f(t,y) (string name, as 'sin')
%         t_interval -- the interval [t0,tn]
%         y0         -- the initial value
%         h          -- the step size
%
% output: t, y       -- the node and the value of y
%
% GP Conangla, 12/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y] = pc_adams4(f_edo, t_interval, y0, h)

t = t_interval(1) : h : t_interval(2);

% use Runge-Kutta method to get four initial values
[~, y] = runge4(f_edo, [t(1),t(4)], y0, h);
y = y';
t = t';

% use predictor - corrector
for n = 4:(length(t)-1)
    % this y(:, n+1) is the predicted value
    y(:, n+1) = y(:,n) + h/24*(55*feval(f_edo, t(n), y(:,n)) - 59*feval(f_edo, t(n-1), y(:,n-1)) + ...
            37*feval(f_edo, t(n-2), y(:,n-2)) - 9*feval(f_edo, t(n-3), y(:,n-3))); % predictor
        
    % now the value y(:, n+1) is updated with the corrected value
    y(:, n+1) = y(:,n) + h/24*(9*feval(f_edo, t(n+1), y(:,n+1)) + 19*feval(f_edo, t(n), y(:,n)) - ...
            5*feval(f_edo, t(n-1), y(:,n-1)) + feval(f_edo, t(n-2), y(:,n-2)));    % corrector
end

% return column vectors with results
t = t';
y = y';
end
