function [iY, qY, y, t] = generateWaveform(pulse, alpha, Ts, n, fc, iData, qData)
%generateWaveform Simulates a received signal.
%   pulse - Function name for pulse generating function
%   alpha - Pulse parameter
%   data - Vector with symbols
%

width = max(length(iData), length(qData));
t = -Ts*width/2:Ts/n:Ts*width/2+(length(iData)-1)*Ts;
iY = pulsetr(pulse, alpha, Ts, n, width, iData);
qY = pulsetr(pulse, alpha, Ts, n, width, qData);
y = iY .* cos(2*pi*fc*t) + qY .* sin(2*pi*fc*t);

% TODO: Add random delay
% TODO: Add noise

end

