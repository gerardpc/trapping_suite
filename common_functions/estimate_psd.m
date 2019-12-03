%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate_psd
% 
% Use Bartlett's method to estimate PSD: divide trace in n_avg segments, 
% calculate the periodogram of each of them and average.
% 
% COMMENTS:
% Output is double sided scaling, but ww \in 2*pi*[0 f_s/2]
% To get single sided scaling, multiply output Sx by 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUTS: 
%               tt      - time vector
%               xx      - vector of trace (i.e., xx(tt))
%               n_avg   - number of averages
%
% OUTPUTS:  
%               Sx      - psd vector
%               ww      - omega vector
%
% Example call:
% >> [Sx, ww] = estimate_psd(tt, xx, 50)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 10.18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Sx, ww] = estimate_psd(tt, xx, n_avg)

%% parameters
n_fft = 2^17;
total_samples = length(tt);
fs = 1/(tt(2) - tt(1));
num_samples = floor(total_samples/n_avg);
% preallocate space
Sx_current = zeros(1, n_fft/2 + 1);

%% Estimate PSD with Bartlett's method
for i = 1:n_avg
%     %% Calculate periodogram (only uses fft function)
%     signal_fft = fft((xx((i - 1)*num_samples + 1:i*num_samples)), n_fft); % the fft itself
%     X2_fft = abs(signal_fft).^2; %signal_fft.*conj(signal_fft);
%     L = length(X2_fft);
%     % proper scaling FFT to estimate PSD, same as periodogram
%     [X_fft,ff] = periodogram(xx((i - 1)*num_samples + 1:i*num_samples), hamming(num_samples), n_fft,fs,'onesided');
%     Sx_current = Sx_current + X2_fft(1:L/2+1)/(fs*num_samples); 
%     ww = 2*pi*fs*(0:L/2)/L;    
%     Sx_current = Sx_current + (X_fft)';
    [X_fft, ff] = periodogram(xx((i - 1)*num_samples + 1:i*num_samples), hamming(num_samples), n_fft, fs,'onesided');    
    Sx_current = Sx_current + (X_fft)'/2;
    ww = ff*2*pi;
end
% average Sx_current to get PSD
Sx = Sx_current/(n_avg);

end