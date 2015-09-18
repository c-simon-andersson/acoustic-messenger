function [iWaveform, qWaveform]  = shift2baseband(waveform, t, fc)
%shift2baseband Shifts the waveform to baseband.
%   waveform - The waveform to shift down.
%   t - Time values for the waveform samples.
%   fc - The carrier frequency.
%   Shifts down the 2-dimensional waveform into its orthogonal sin and cos
%   baseband signals.

iWaveform = waveform * sqrt(2) .* cos(2*pi*fc*t);
qWaveform = waveform * sqrt(2) .* sin(2*pi*fc*t);

end

