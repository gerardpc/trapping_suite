function[] = test_sde(n_samples)
n = 2;        % Order of system
t0 = 0;       % Initial time
dt = 1e-3;    % Fixed time step
tf = 1;       % Final time
t = t0:dt:tf; % Time vector
t_length = length(t);

A = 1;
B = 1e-1;
C = 0;

f = @(t,x)[x(2); -(A*x(2) + C*x(1))]; % Drift function
ep = B;                                % Size of additive noise
g = @(t,x)[0;ep];                         % Diffusion function


x = zeros(t_length,n); % Allocate output
var_x = zeros(t_length,n); % Allocate output
x0 = [0;0];
x(1,:) = x0;     % Set initial condition

% %seed = 1;  % Seed value
% rng(seed); % Always seed random number generator

% Euler-Maruyama
for j = 1:n_samples;
    fprintf('%d\n', j);
    for i = 1:t_length-1
        x(i+1,:) = x(i,:).' + f(t(i),x(i,:))*dt + g(t(i),x(i,:)).*randn(n,1)*sqrt(dt);
    end
    var_x = var_x + x.^2;
end
var_x = var_x/(n_samples - 1);
% check
var_short = B^2/(2*A)*(1 - exp(-2*A*t));

figure(1);
clf;
hold on;
plot(t, var_x(:,2));
plot(t, var_short);
xlabel('Time');
ylabel('State');
legend('x(t)','xdot(t)')
end