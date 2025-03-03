function [freq, amp] = get_spec(Fs, data)
% Fs, sample rate in Hz
% data vector
% modified from Jana's code

M=length(data);
NFFT = 2^nextpow2(M); % Next power of 2 from length of data
fO = (Fs)/2*linspace(0,1,NFFT/2+1); % freq. bins, 0 to Nyquist
func=data-mean(data);   % remove DC
BF = fft(func,NFFT);    % FFT
BF=(abs(BF(1:NFFT/2+1)));   % amplitude of FFT;use only first half due to symmetry
freq=fO(1:end);
amp=abs(BF(1:end))./Fs;