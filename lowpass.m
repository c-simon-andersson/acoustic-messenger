function [iWaveformFiltered, qWaveformFiltered] = lowpass(iWaveform, qWaveform, n)
%lowpass Lowpass filter the signal after bringing it down to baseband.
%   iWaveform - The first waveform.
%   qWaveform - The second waveform.
%   n - Samples per symbol.
%   Passes the signals through a MATLAB-generated low-pass filter.

d = fdesign.lowpass('Fp,Fst,Ap,Ast', 1/n,2/n,0.1,60);
Hd = design(d, 'butter');
iWaveformFiltered = filter(Hd, iWaveform);
qWaveformFiltered = filter(Hd, qWaveform);

end

