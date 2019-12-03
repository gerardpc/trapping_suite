%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full_numerical_ode.m
%
% Numerically solve ODE with Adams 4th order predictor corrector method
% ODE type: y' = f(t,y), where f is saved in function test.m/test_2nd.order.m
% Output is plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS: 
%       time_interval  -- solve y(t) for t \in [a b]
%       h              -- discretization size
%       order          -- 2nd/1st order ode (eg 'order', '1')
%       initial_values'-- initial values (eg 'ini_values', '[1 0]')
%       'PSD' := Calculate or not PSD of ode solution ('PSD', 'on'/'off')
% OUTPUTS:
%       t              -- Time vector
%       y              -- vector with numerical solution (if 1st order, 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GP Conangla 12/2015, updated 7/2018 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[t, y] = full_numerical_ode(time_interval, h, order, initial_values)

%% variable arguments in
switch(order)
    case '1st_order'
        order = 'test_1st_order';
    case '2nd_order'
        order = 'test_2nd_order';
end

%% Solve numerically
[t, y] = pc_adams4(order, time_interval, initial_values, h);     

%% Calculate plot

% Plot result
figure(1);
clf;
subplot(2,1,1);
nice_plot(t, y(:,1), 'time', 'x(t)', 'ODE solution by 4th order PC Adams method', 'color', 'r');
subplot(2,1,2);
nice_plot(t, y(:,2), 'time', 'v(t)', 'ODE solution by 4th order PC Adams method');

end