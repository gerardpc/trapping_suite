

%% Plot results (PSD + fit)
figure(1);
clf;
[Sx, ww] = estimate_psd(t, x_kf(1,:), 100);
nice_plot(ww/2/pi, Sx, '$f$ (Hz)', '$S_x$', 'Estimated PSD', 'semilogy', 'linewidth', '1');
area(ww/2/pi, Sx); 
[Sx, ww] = estimate_psd(t, x(1,:), 100);
nice_plot(ww/2/pi, Sx, '$f$ (Hz)', '$S_x$', 'Estimated PSD', 'semilogy', 'linewidth', '1');
area(ww/2/pi, Sx); 
[Sx, ww] = estimate_psd(t, x_1e3_nf(1,:), 100);
nice_plot(ww/2/pi, Sx, '$f$ (Hz)', '$S_x$', 'Estimated PSD', 'semilogy', 'linewidth', '1');
area(ww/2/pi, Sx); 
xlim([0 5e5]);
ylim([1e-29 1e-16]);