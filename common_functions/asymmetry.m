[time,x_mode] = extract_trace('01_TimeTrace_delay=210_kd=-100.bin', '\\nas01\GroupPrivate\PNO\MEMBERS\Francesco\Data\ElectricCooling\P15_235nm\PSD_Asymmetry', 20);

[Sx, ww] = estimate_psd(time, x_mode(4,:), 1000);
figure(1);
clf;
nice_plot(ww/2/pi/1e3, Sx, 'freq (kHz)', 'S_x(\omega)', 'PSD asymmetry', 'semilogy');
xlim([105 145]);