function [p,P,f] = applycalcurve(fs,v,f_cc,F_cc)
% This function takes a voltage waveform (in time domain) and multiplies it
% by frequency-dependent calibration factors in frequency domain, returning
% pressure values. Inter/Extrapolation of the calibration curve is used, so
% make sure the curve is defined throughout the spectrum of the waveform.
%
% Input:
%   - fs: sampling rate (Hz)
%   - v: input signal (V)
%   - f_cc: frequencies at which calibration curve is defined (Hz)
%   - F_cc: calibration factors as a function of frequency (Pa/Hz)
%
% Output:
%   - p: pressure signal, in time domain (Pa)
%   - P: pressure signal, in frequency domain (Pa)
%   - f: frequencies at which the spectrum is defined (Hz)
%
%%% Reynard de Vries
%%% TU Delft
%%% First version: 05 May 2016
%%% Last update: 20 June 2016

%% Fourier transform input
% Length of input signal 
NFFT = length(v);  

% Fourier coefficients (V)
V_2 = fft(v,NFFT)/NFFT;     

% Define frequency vector (Hz)
f = 0:fs/NFFT:fs/2;     

% One-sided Fourier coefficients (Pa)
V = V_2(1:length(f));       


%% Multiply voltage spectrum by calibration curve
    
% Interpolate cal curve at specified frequencies (Real number)
F_int = interp1(f_cc,F_cc,f,'linear','extrap');

% Multiply point-by-point (Complex number times real number)
P = transpose(V).*F_int;

%% Inverse Fourier transform
% Two-sided Fourier pressure coefficients (Pa)
P_2     = [P flipud(conj(P(2:ceil(NFFT/2))))]; 

% Corrected signal data in time domain (Pa)
p = ifft(P_2,NFFT,'symmetric')*NFFT; 


%% Plots for checking
% figure
% semilogx(f,V,'b');
% hold on;
% semilogx(f,P,'r');
% legend('V(f)','P(f)');
% xlabel('Frequency')
% grid on;
% 
% figure
% plot(v(1:500),'b');
% hold on;
% plot(p(1:500),'r');
% legend('v(t)','p(t)');
% xlabel('Time')
% grid on;
    

