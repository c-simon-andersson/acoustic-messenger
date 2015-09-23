function [iWaveform, qWaveform]  = shift2baseband(waveform, t, fc, Ts, n)
%shift2baseband Shifts the waveform to baseband.
%   waveform - The waveform to shift down.
%   t - Time values for the waveform samples.
%   fc - The carrier frequency.
%   Shifts down the 2-dimensional waveform into its orthogonal sin and cos
%   baseband signals.

% Shift to baseband using uncertain carrier frequency
iWaveform = waveform * sqrt(2) .* cos(2*pi*fc*t);
qWaveform = waveform * sqrt(2) .* sin(2*pi*fc*t);

% Pass the baseband signal through a low-pass filter
% TODO: Filter should be created on startup instead of every time we want
% to use it. It's computationally expensive to create the filter.
[iWaveform, qWaveform] = lowpass(iWaveform, qWaveform, n);

% Find second peak using early late (skip first to avoid transients)
epsilon = 0.000001 * max(iWaveform);
epsilon = [epsilon epsilon];
ns = 10;

% Find the middle chunk of the square wave
% TODO: Improve the detection of the beginning and end of the square wave.
% Current implementation has problems with short pilot tones.
c = earlyLate(iWaveform, ns, 1, epsilon);
d = earlyLate(flip(iWaveform), ns, 1, epsilon);
d = length(iWaveform) - d;
middlePoint = floor(c + (d - c) / 2);
startPoint = floor(middlePoint - (d - c) / 10);
endPoint = floor(middlePoint + (d - c) / 10); 

% Sum up the angles between consecutive points
% TODO: This method is computationally expensive. Can we find some more
% efficient solution?
angle = 0;
for currentPoint = startPoint:(endPoint - 1)
    v1 = [iWaveform(currentPoint) qWaveform(currentPoint)];
    v2 = [iWaveform(currentPoint + 1) qWaveform(currentPoint + 1)];
    angle = angle + acos((v1 * v2') / (sqrt(sum(v1.^2)) * sqrt(sum(v2.^2))));   
end

% Calculate frequency
angle = real(angle);
time = (endPoint-startPoint)/n*Ts;
period = time / (angle/(2*pi));
frequency = 1/period;

% Calculate sign
% TODO: Doing this over just one sample is sensitive to noise. Should be
% done over the entire range instead.
delta = iWaveform(startPoint)*qWaveform(startPoint + 1) - iWaveform(startPoint + 1)*qWaveform(startPoint);
if delta > 0
   frequency = frequency * -1;
end


figure
grid on
hold on
plot(t, iWaveform)
plot(t, qWaveform)
plot(t(middlePoint), iWaveform(middlePoint), 'ko')
plot(t(startPoint), iWaveform(startPoint), 'bo')
plot(t(endPoint), iWaveform(endPoint), 'bo')

% Revise carrier frequency
fc = fc + frequency;

% Shift to baseband using calculated carrier frequency
iWaveform = waveform * sqrt(2) .* cos(2*pi*fc*t);
qWaveform = waveform * sqrt(2) .* sin(2*pi*fc*t);

% Pass the baseband signal through a low-pass filter
[iWaveform, qWaveform] = lowpass(iWaveform, qWaveform, n);

% Find the i and q values
iAvg = sum(iWaveform(startPoint:endPoint)) / (endPoint - startPoint);
qAvg = sum(qWaveform(startPoint:endPoint)) / (endPoint - startPoint);

% Calculate phase
v1 = [iAvg qAvg];
v2 = [1 1]; % Reference vector
phase = acos((v1 * v2') / (norm(v1) * norm(v2)))
delta = v1(1)*v2(2) - v1(2)*v2(1);
if delta > 0
   phase = phase * -1; 
end

% Shift to baseband using calculated carrier frequency
iWaveform = waveform * sqrt(2) .* cos(2*pi*fc*t - phase);
qWaveform = waveform * sqrt(2) .* sin(2*pi*fc*t - phase);

% Pass the baseband signal through a low-pass filter
[iWaveform, qWaveform] = lowpass(iWaveform, qWaveform, n);

plot(t, iWaveform, 'LineWidth', 2)
plot(t, qWaveform, 'LineWidth', 2)

end

