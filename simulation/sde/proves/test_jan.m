
function[] = test_jan(t_end, h, n_averages)
% Definitions for each simulation:
% t_end: run simulation until t = t_end
% n_averages = number of averages to estimate moments (later will be used to apply Richardson extrapolation)

% Definitions: experimental parameters (like in pdf)
r = 117.5e-9; % radius particle
density = 1850; % kg/m^3, silica
m = 4/3*pi*r^3*density;
gamma_air = 6*pi*r*1.84e-5; % stokes law
pressure = 1000; % in mbar
gamma = gamma_air*pressure/1010;
Gamma = gamma/m;
eps = 6.3e-9;
T = 295;
k_B = 1.38065e-23;
sigma = sqrt(2*k_B*T*gamma);
w_0 = 2*pi*150e3;

% % Definitions: experimental parameters (like in pdf)
% m = 9.2e-18;
% gamma = 3.5e-11;
% w = 2*pi*2e4;
% eps = 6.3e-9;
% T = 295;
% k_B = 1.38065e-23;
% sigma = sqrt(2*k_B*T*gamma);

% % Parameters as defined in Jan & Dwight's calculations
% lambda_1 = -m*eps^2/(2*gamma^3);
% tau = 1/lambda_1;

% Run simulation
[t, y, mean_y, var_y] = full_numerical_sde([0 t_end], h, '2nd_order', [0 0], n_averages);

% % Compare with theory
% t = linspace(0, t_end, 1e3);
% w_1 = sqrt(w_0^2 - Gamma^2/4);
% var_x_opt = k_B*T/(m*w_0^2)*(1 - exp(-Gamma*t).*(w_0^2/w_1^2 - Gamma^2/(2*w_1)^2*cos(2*w_1*t) + Gamma/(2*w_1)*sin(2*w_1*t)));
% var_eq = k_B*T/(m*w_0^2);
% 
% % Plot results
% % plot var estimation
% figure(2);
% clf;
% box;
% hold on;
% % Add simulation result to plot
% plot(t_2, var_y_2(:,1), 'Color', 'b', 'LineWidth', 3.5);
% plot(t, ones(1,length(t))*var_eq, '--', 'Color', 'r', 'LineWidth', 3.5);
% plot(t, var_x_opt, 'Color', 'r', 'LineWidth', 3.5);
% set(gca,'FontSize',18,'FontName', 'CMU Sans Serif');
% xlabel('Time', 'fontsize', 24, 'FontName', 'CMU Sans Serif');
% ylabel('$\bf{E}[y^2]$', 'fontsize', 24, 'Interpreter', 'latex', 'FontName', 'CMU Sans Serif');
% title('Euler-Maruyama method');
% xlim([t(1), t(end)]);
% grid on;

% Compare with theory
A = 1;
B = 1e-1;
C = 2;
C2 = (C - A^2/4);
%var_short = B^2/(2*A)*(1 - exp(-2*A*t_2));
var_short = B^2/(2*A*2)*(1 - exp(-A*t).*(C/C2 - A^2/(4*C2)*cos(2*sqrt(C2)*t) + A/(2*sqrt(C2))*sin(2*sqrt(C2)*t)));
% Plot results
% plot var estimation
figure(2);
clf;
box;
hold on;
% Add simulation result to plot
plot(t, var_y(:,1), 'Color', 'b', 'LineWidth', 3.5);
plot(t, var_short, '--', 'Color', 'r', 'LineWidth', 3.5);
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif');
xlabel('Time', 'fontsize', 24, 'FontName', 'CMU Sans Serif');
ylabel('Estimated $\bf{E}[y^2]$', 'fontsize', 24, 'Interpreter', 'latex', 'FontName', 'CMU Sans Serif');
string_title = sprintf('%d runs of Runge-Kutta, 1st order', n_averages);
title(string_title);
xlim([t(1), t(end)]);
grid on;
