function [pack, psd, const, eyed] = receiver_v2(tout,fc)
% Receiver with barker sequence indentification and phase correction

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

a = 0.35; tau = 1/sym_rate; span = 4; %a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
rrc_pulse = rtrcpuls(a,tau,fs,span);
match_filter = rrc_pulse;

%% Filter construction
LP_filter=firpm(30, [0 0.50 0.55 1], [1 1 0 0]); % TODO: Tune filter parameters
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