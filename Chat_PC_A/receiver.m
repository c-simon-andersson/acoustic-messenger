function [pack, psd, const, eyed] = receiver(tout,fc)
%receiver PC_A receiver script
%   tout - Time out for receiver. Will return after this time if no message
%   is detected.
%   fc - Carrier frequency to listen on.

pack = []; psd = [];  const=[]; eyed = [];

recTime = 1; % Seconds to record waveform
alpha = 0.3; % Roll-off
pulse = 'rtrcpuls';
pulse = 'rcpuls';
Ts = 1/240; % Symbol time

% Data to transmit. Used to generate waveform during receiver development.
qData = [1 0 0 1 0 0];
iData = [0 1 0 0 0 0];

f2 = (1+alpha)/(2*Ts); % Pulse bandwidth
n = max(4*(fc+f2/2),4*f2); % Samples per symbol

% Wait for a message for tout amount of time.
foundMessage = detectMessage(tout, fc);

% Exit if we didn't find anything.
if foundMessage == 0
   return;
end

% Record the raw waveform
%noisyWaveform = recordWaveform(recTime);
[iY, qY, noisyWaveform, t] = generateWaveform(pulse, alpha, Ts, n, fc, iData, qData);

figure
hold on
grid on
plot(t, iY, 'LineWidth', 2)
plot(t, qY, 'LineWidth', 2)
plot(t, noisyWaveform, '--')

% Frequency/phase detection and shift to baseband
[iWaveform, qWaveform] = shift2baseband(noisyWaveform, t, fc);

% Pass the baseband signal through a low-pass filter
[iWaveformFiltered, qWaveformFiltered] = lowpass(iWaveform, qWaveform, n);

plot(t, iWaveformFiltered)
plot(t, qWaveformFiltered)

% Run waveform through the matched filter
[iMatched, qMatched] = matchedFilter(iWaveformFiltered, qWaveformFiltered, pulsetr(pulse, alpha, Ts, n, 2,1));

figure
hold on
grid on
t = 1:length(iMatched);
plot(t, iMatched)
plot(t, qMatched)

% Detect sampling time and sample
symbolSequence = sampleMatched(iMatched, qMatched);

% Frame synchronization and convert symbol sequence to bits
bitSequence = symbols2bits(symbolSequence);

pack = bitSequence;


end