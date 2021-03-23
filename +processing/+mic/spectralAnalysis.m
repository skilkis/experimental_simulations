%% Function spectralAnalysis.m
% Converts input time domain signal to the frequency domain using Welch's 
% method 
% =========================================================================
% Tomas Sinnige - T.Sinnige@tudelft.nl 
% TU Delft - LR - AWEP - Flight Performance and Propulsion
%
% Version: 1.0
% Last updated:  22 April 2015
% First version: 22 April 2015
% =========================================================================
% | Version |    Date   |   Author  |              Changelog              |
% |---------|-----------|-----------|-------------------------------------|
% |   1.0   | 22/04/'15 | T.Sinnige | First version                       |
% |---------|-----------|-----------|-------------------------------------|
% =========================================================================
% Inputs: y      - signal to be converted to frequency domain. 
%                    !!time axis should be along first dimension !!
%         w      - window function
%         wOvrlp - overlap window function
%         fS     - sampling frequency of input signal y [Hz]
% -------------------------------------------------------------------------  
% Output: PXX - PSD of input signal [unit/Hz]
%         SPL - sound pressure levels computed from PXX [dB]
%         f   - frequency vector [Hz] at which the ydB values are defined
% =========================================================================
function [PXX,SPL,f] = spectralAnalysis(y,w,wOvrlp,fS,fEval)

%% Compute window properties
N  = length(w);                  % number of elements per block
CG = (1/N)*sum(w(1:end-1,:));    % coherent gain
NG = (1/N)*sum(w(1:end-1,:).^2); % noise gain

%% Use pwelch to convert into frequency domain
if verLessThan('matlab','8.3') % if Matlab version older than 2014a, the pwelch function has to be called in a loop
    for i=1:size(y,2)
        [pxx(:,i),f] = pwelch(y(:,i),w,wOvrlp,fEval,fS);
    end
else % if newer Matlab version, the pwelch function call can be done at once for the whole matrix
    [pxx,f] = pwelch(y,w,wOvrlp,fEval,fS);
end

%% Compute power spectral density
fbin = f(2)-f(1); % frequency bin width [Hz]

% PXX = pxx*fbin*NG/CG^2; % power spectral density [unit/Hz] (pxx multiplied with frequency bin width to get power of signal, see http://www.fhnw.ch/technik/ime/publikationen/2012/how-to-use-the-fft-and-matlab2019s-pwelch-function-for-signal-and-noise-simulations-and-measurements)
% display('coherent gain correction applied --> proper tonal amplitude results')

% PXX = pxx; % power spectral density [unit/Hz] (pxx multiplied with frequency bin width to get power of signal, see http://www.fhnw.ch/technik/ime/publikationen/2012/how-to-use-the-fft-and-matlab2019s-pwelch-function-for-signal-and-noise-simulations-and-measurements)
% display('coherent gain correction NOT applied --> proper background noise')

PXX = pxx*fbin; % power spectral density [unit/Hz] (pxx multiplied with frequency bin width to get power of signal, see http://www.fhnw.ch/technik/ime/publikationen/2012/how-to-use-the-fft-and-matlab2019s-pwelch-function-for-signal-and-noise-simulations-and-measurements)
display('frequency bin width correction applied --> proper spectral levels for integration')

%% Compute Sound Pressure Level from PSD
p0  = 20e-6; % reference sound pressure [Pa] 
SPL = 10*log10(PXX/p0^2); % sound pressure level [dB]
end % end of function APIANINF_spectralAnalysis.m