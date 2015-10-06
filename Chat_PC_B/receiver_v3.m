function [pack, psd, const, eyed] = receiver_v3(tout,fc)
% Receiver with barker sequence indentification and phase correction

pack = []; psd = [];  const=[]; eyed = [];

%% Definitions
barker = [1 1 1 0 0 0 1 0 0 1 0];
pilot = ones(1, 30);

wave = load('wave.mat');
wave = wave.output;

n_bits = 416;
syms_per_bit = 2;
sym_rate = 240;
fs = 48e3;
t = (1:numel(wave))/fs;
pl = floor(length(pilot)*3/4) %pilot length in symbols

a = 0.35; tau = 1/sym_rate; span = 4; %a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
rrc_pulse = rtrcpuls(a,tau,fs,span);
match_filter = rrc_pulse;

%% Filter construction
LP_filter=firpm(30, [0 0.50 0.55 1], [1 1 0 0]); % TODO: Tune filter parameters
d = fdesign.lowpass('Fp,Fst,Ap,Ast', 1/(fs/sym_rate),2/(fs/sym_rate),0.1,60);
Hd = design(d, 'butter');
LPMF = conv(match_filter, LP_filter); % Combining LP and match-filter

%% Synchronization
% Construct barker sequence filter
barker_upsampled = upsample(barker*2-1, fs/sym_rate);
barker_filter = conv(barker_upsampled, rrc_pulse);

barker_center = zeros(1,24);
barker_val = zeros(1,24);
for i = 1:24;
    phase_shift = i * 2*pi/24;
    MFout_real = conv(wave.*cos(2*pi*fc*t + phase_shift)*sqrt(2), LPMF); 
    MFout_imag = conv(wave.*sin(2*pi*fc*t + phase_shift)*sqrt(2), LPMF);

    % Convolve the signals to find maximum correlation.
    barker_signal_real = fliplr(conv(fliplr(MFout_real), barker_filter, 'same'));
    barker_signal_imag = fliplr(conv(fliplr(MFout_imag), barker_filter, 'same'));

    barker_signal_sum = barker_signal_real+barker_signal_imag;
    [barker_val(i),barker_center(i)] = max(barker_signal_sum);
end
[~,phase_ind] = max(barker_val);
signal_start = barker_center(phase_ind) + (numel(barker)/2)*fs/sym_rate;

phase_shift = phase_ind * 2*pi/24;
MFout_real = conv(wave.*cos(2*pi*fc*t + phase_shift)*sqrt(2), LPMF); 
MFout_imag = conv(wave.*sin(2*pi*fc*t + phase_shift)*sqrt(2), LPMF);

MFout_real = conv(filter(Hd, wave.*cos(2*pi*fc*t + phase_shift)*sqrt(2)), match_filter); 
MFout_imag = conv(filter(Hd, wave.*sin(2*pi*fc*t + phase_shift)*sqrt(2)), match_filter);
signal_start = signal_start+824; %filter delay
%[MFout_real, MFout_imag] = lowpass(MFout_real, MFout_imag, fs/sym_rate);
%[MFout_real, MFout_imag] = shift2baseband(wave, t, fc, 1/sym_rate, fs/sym_rate, pl*2, signal_start);

currentPoint = signal_start;
% Set the start and end point to calculate angle on
angleInterval = floor(fs/sym_rate*pl); % Length of interval to measure angle on
startPoint = currentPoint; endPoint = startPoint + angleInterval;

% Sum up the angles between consecutive points
% TODO: Probably would be nicer to do this using vectors
angle = 0;
for currentPoint = startPoint:(endPoint - 1)
    v1 = [MFout_real(currentPoint) MFout_imag(currentPoint)];
    v2 = [MFout_real(currentPoint + 1) MFout_imag(currentPoint + 1)];
    thisAngle = acos((v1 * v2') / (sqrt(sum(v1.^2)) * sqrt(sum(v2.^2))));
    delta = v1(1)*v2(2)-v1(2)*v2(1);
    if delta > 0
       thisAngle = thisAngle * -1; 
    end
    angle = angle + thisAngle;
end

% Calculate frequency
angle = real(angle);
time = (endPoint-startPoint)/(fs/sym_rate)*1/sym_rate;
period = time / (angle/(2*pi));
frequency = 1/period

figure
grid on
hold on
plot(MFout_real)
plot(MFout_imag)
plot(startPoint, MFout_real(startPoint), 'bo')
plot(endPoint, MFout_real(endPoint), 'bo')
plot(signal_start, MFout_real(signal_start), 'ko')

fc = fc + frequency;
%MFout_real = conv(wave.*cos(2*pi*fc*t + phase_shift)*sqrt(2), LPMF); 
%MFout_imag = conv(wave.*sin(2*pi*fc*t + phase_shift)*sqrt(2), LPMF);

MFout_real = conv(filter(Hd, wave.*cos(2*pi*fc*t + phase_shift)*sqrt(2)), match_filter); 
MFout_imag = conv(filter(Hd, wave.*sin(2*pi*fc*t + phase_shift)*sqrt(2)), match_filter);
%[MFout_real, MFout_imag] = lowpass(MFout_real, MFout_imag, fs/sym_rate);


%% Sampling
sample_vec = zeros(1, n_bits/syms_per_bit);

sample_vec(1) = signal_start;
for i = 2:numel(sample_vec)
    sample_vec(i) = sample_vec(i-1) + fs/sym_rate;
end
const = [MFout_real(sample_vec); MFout_imag(sample_vec)];

%% Plotting
figure
subplot(2,1,1)
plot(MFout_real)
title('MF output, real')

subplot(2,1,2)
plot(MFout_imag)
title('MF output, imag')
% 
% figure
% subplot(2,1,1)
% plot(barker_filter)
% title('barker filter')
% 
% subplot(2,1,2)
% plot(barker_signal_sum)
% title('barker signal sum')
% hold on
% plot(barker_center, barker_signal_sum(barker_center), 'o')


%% Output
%pack = symbols2bits(const', [1 1 1 0 0 1 0]);
pack = (sign(reshape(const, 1, n_bits))+1)/2;
end