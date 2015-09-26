function [pack, psd, const, eyed] = receiver(tout,fc)
% Receiver with barker sequence indentification
% Cannot handle phase shift - check receiver_v2

pack = []; psd = [];  const=[]; eyed = [];

%% Definitions
barker = [1 1 1 0 0 0 1 0 0 1 0];

wave = load('wave.mat');
wave = wave.output;

n_bits = 416;
syms_per_bit = 2;
sym_rate = 240;
fs = 48e3;
t = (1:numel(wave))/fs;

%a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4;
rrc_pulse = rtrcpuls(a,tau,fs,span);
match_filter = rrc_pulse;

%% Filter construction and filtering
LP_filter=firpm(30, [0 0.50 0.55 1], [1 1 0 0]); % TODO: Tune filter parameters
LPMF = conv(match_filter, LP_filter); % Combining LP and match-filter

MFout_real = conv(wave.*cos(2*pi*fc*t)*sqrt(2), LPMF); 
MFout_imag = conv(wave.*sin(2*pi*fc*t)*sqrt(2), LPMF);

figure
subplot(2,1,1)
plot(MFout_real)
title('MF output, real')

subplot(2,1,2)
plot(MFout_imag)
title('MF output, imag')

%% Synchronization
% Construct barker sequence filter
barker_upsampled = upsample(barker*2-1, fs/sym_rate);
barker_filter = conv(barker_upsampled, rrc_pulse);

% Convolve the signals to find maximum correlation.
barker_signal_real = fliplr(conv(fliplr(MFout_real), barker_filter, 'same'));
barker_signal_imag = fliplr(conv(fliplr(MFout_imag), barker_filter, 'same'));

% % This part has no solid theoretical ground that I know of
% barker_sum = sum(abs(barker_filter));
% barker_signal_1 = barker_signal_real + barker_signal_imag;
% barker_signal_2 = barker_signal_1;
% barker_signal_2(barker_signal_2 < barker_sum/2) = 0;

barker_signal_sum = barker_signal_real+barker_signal_imag;
[~,barker_center] = max(barker_signal_sum);
% find(diff(barker_signal_2)<0, 1, 'first');

signal_start = barker_center + (numel(barker)/2)*fs/sym_rate;

figure
subplot(2,1,1)
plot(barker_filter)
title('barker filter')

subplot(2,1,2)
plot(barker_signal_sum)
title('barker signal sum')
hold on
plot(barker_center, barker_signal_sum(barker_center), 'o')

%% Sampling
sample_vec = zeros(1, n_bits/syms_per_bit);

sample_vec(1) = signal_start;
for i = 2:numel(sample_vec)
    sample_vec(i) = sample_vec(i-1) + fs/sym_rate;
end
const = [MFout_real(sample_vec); MFout_imag(sample_vec)];

figure

subplot(2,1,1)
hold on
plot(MFout_real)
plot(sample_vec, MFout_real(sample_vec),'o')
title('MF output, real')

subplot(2,1,2)
hold on
plot(MFout_imag)
plot(sample_vec, MFout_imag(sample_vec),'o')
title('MF output, imag')


%% Output
%pack = symbols2bits(const', [1 1 1 0 0 1 0]);
pack = (sign(reshape(const, 1, n_bits))+1)/2;
end