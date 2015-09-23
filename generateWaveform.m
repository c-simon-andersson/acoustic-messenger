function [iY, qY, y, t] = generateWaveform(pulse, alpha, Ts, n, fc, phase, iData, qData)
%generateWaveform Simulates a received signal.
%   pulse - Function name for pulse generating function
%   alpha - Pulse parameter
%   data - Vector with symbols
%

width = 4*max(length(iData), length(qData));
t = -Ts*width/2:Ts/n:Ts*width/2+(length(iData)-1)*Ts;
iY = pulsetr(pulse, alpha, Ts, n, width, iData);
qY = pulsetr(pulse, alpha, Ts, n, width, qData);
y = iY .* cos(2*pi*fc*t + phase) + qY .* sin(2*pi*fc*t + phase);

% TODO: Add random delay

% Add AWGN
%y = y + awgn(y, 100);


end

