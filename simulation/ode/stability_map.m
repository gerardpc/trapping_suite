
x = 0:0.2:3;
y = 0:0.2:2;
[Q, Gamma] = meshgrid(x, y);
Stab_map = 0.01*ones(size(Q));
global q;
global gamma;
for i = 1:length(y)
    for j = 1:length(x);
        q = x(j);
        gamma = y(i);
        Stab_map(i, j) = full_numerical_ode([0 100], 1e-2, '2nd_order', 'PSD_off');
        fprintf('q = %.2g,  gamma = %.2g\n', q, gamma);
        figure(1);
        surf(Q, Gamma, log(Stab_map));
        %plot(Q, (Stab_map));
        az = 0;
        el = 90;
        view(az, el);        
        drawnow;
    end
end
