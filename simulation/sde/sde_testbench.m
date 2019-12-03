%% Test Bench: Compare with theory
% Initial conditions = 0
%
% Some parameters: 
% A = 1;
% B = 1e-1;
% C = 2;
%
% O-U, 1st order
% f(1) = -A*x; 
% sigma = B
% Then, variance of velocity will be
% var_short = B^2/(2*A)*(1 - exp(-2*A*t));
%
% Free Brownian motion, 2nd order
% f(1) = v; 
% f(2) = -A*v - C*x;
% sigma = B
% Then, variance of velocity will be
% var_short = B^2/(2*A)*(1 - exp(-2*A*t));
%
% Stochastic harmonic oscillator, 2nd order% f(1) = v; 
% f(2) = -A*v - C*x;
% sigma = B
% Then, variance of x will be
% var_short = B^2/(2*A*C)*(1 - exp(-A*t).*(C/C2 - A^2/(4*C2)*cos(2*sqrt(C2)*t) + A/(2*sqrt(C2))*sin(2*sqrt(C2)*t)));
% where C2 = (C - A^2/4);

%% Run simulation
% t_end = 5;
% h = 1e-2;
% n_averages = 2000;
%
% [t, y, mean_y, var_y] = full_numerical_sde([0 t_end], h, order, [0 0], n_averages);
%
%% Plot results of testbench
% figure;
% clf;
% box;
% hold on;
% % Add simulation result to plot
% plot(t, var_y(:,1), 'Color', 'b', 'LineWidth', 3.5);
% plot(t, var_short, '--', 'Color', 'r', 'LineWidth', 3.5);
% set(gca,'FontSize',18,'FontName', 'CMU Sans Serif');
% xlabel('Time', 'fontsize', 24, 'FontName', 'CMU Sans Serif');
% ylabel('Estimated $\bf{E}[y^2]$', 'fontsize', 24, 'Interpreter', 'latex', 'FontName', 'CMU Sans Serif');
% string_title = sprintf('%d runs of Runge-Kutta, 1st order', n_averages);
% title(string_title);
% xlim([t(1), t(end)]);
% grid on;