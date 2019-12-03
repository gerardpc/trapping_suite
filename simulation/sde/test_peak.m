
w_dr = 2*pi*130e3;
f_0 = 10*1.60217662e-19*3.7/2*564;
tau = 20e-3;
opt = load_opt_tweezer();
m = opt.m;
w_0 = opt.w_0;
G_norm = opt.G_norm;

peak_height = f_0^2*tau/(4*m^2*((w_dr^2 - w_0^2)^2 + G_norm^2*w_dr^2));
fprintf('Height of peak (analytical): %.3e\n', peak_height);

% % vec_peak = zeros(10, 1);
% % for i = 1:10
% %     vec_peak(i) = max(mean(psdcell{i}, 2));
% % end
% % vec_peak = vec_peak - opt.sigma^2/(m^2*((w_dr^2 - w_0^2)^2 + G_norm^2*w_dr^2));
% % fprintf('Height of peak from simulation: %.3e\n', mean(vec_peak));
% % fprintf('Std of peak from simulation:    %.3e\n', std(vec_peak));
% 
% for i = 1:6
%     for j = 1:5
%         position = i + j*10;
%         smooth_psd = mean(psdcell_2{position}, 2);
%         psdcell_2(position) = [];
%         psdcell_3{position} = smooth_psd;
% %         figure(position);
% %         clf;
% %         nice_plot(f_period, (smooth_psd), '$f$ (Hz)', '$S_x$', 'Estimated PSD from Monte Carlo traces', 'semilogy', 'linewidth', '1');
% %        xlim([0 2*opt.w_0/(2*pi)]);
%         pause(1);
%     end
% end