function [pack, psd, const, eyed] = receiver(tout,fc)
% Receiver with barker sequence indentification and phase correction

%%%% Definitions
[pack, psd, const, eyed] = deal([]);
barker = [1 1 1 0 0 0 1 0 0 1 0];
n_bits = 432;
syms_per_bit = 2;
sym_rate = 120;
fs = 12e3;
rec_bits = 16;

barker_threshold = 120;
phase_resolution = 12;
wave_start = 1;

% a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4; 
rrc_pulse = rtrcpuls(a,tau,fs,span);
match_filter = rrc_pulse;

barker_upsampled = upsample(barker*2-1, fs/sym_rate);
barker_filter = conv(barker_upsampled, rrc_pulse);

%%%% Filter construction
% LP_filter = firpm(30, [0 0.50 0.55 1], [1 1 0 0]); % TODO: Tune filter parameters
% LPMF = conv(match_filter, LP_filter); % Combining LP and match-filter
LPMF = match_filter;


rec = audiorecorder(fs, rec_bits, 1);
record(rec);
pause(0.1)
tic;
%%%% main loop
disp('recieving')
while toc < tout && isempty(pack)

    wave = getaudiodata(rec, 'int16');
    wave_end = numel(wave);
     
    wave = double(wave(wave_start:end)');
    wave = wave/max(wave);
    t = (1:numel(wave))/fs;
    if( numel(wave) < 1.1*numel(barker_filter) )
        continue;
    end
    
    %%%% Synchronization
    barker_center = zeros(1,phase_resolution);
    barker_val = zeros(1,phase_resolution);
    for i = 1:phase_resolution;
        phase_shift = i * 2*pi/phase_resolution;
        MFout_real = conv(wave.*cos(2*pi*fc*t + phase_shift)*sqrt(2), LPMF); 
        MFout_imag = conv(wave.*sin(2*pi*fc*t + phase_shift)*sqrt(2), LPMF);

        % Convolve the signals to find maximum correlation.
        barker_signal_real = fliplr(conv(fliplr(MFout_real), barker_filter, 'same'));
        barker_signal_imag = fliplr(conv(fliplr(MFout_imag), barker_filter, 'same'));

        barker_signal_sum = barker_signal_real+barker_signal_imag;
        [barker_val(i),barker_center(i)] = max(barker_signal_sum);
    end
    [max_barker,phase_ind] = max(barker_val);
    
    if max_barker < barker_threshold && ~exist('sample_vec','var')
        wave_start = wave_end;
        continue;
    end

    signal_start = barker_center(phase_ind) + (numel(barker)/2)*fs/sym_rate;

    phase_shift = phase_ind * 2*pi/phase_resolution;
    MFout_real = conv(wave.*cos(2*pi*fc*t + phase_shift)*sqrt(2), LPMF); 
    MFout_imag = conv(wave.*sin(2*pi*fc*t + phase_shift)*sqrt(2), LPMF);


    %%%% Sampling
    sample_vec = zeros(1, n_bits/syms_per_bit);

    sample_vec(1) = signal_start;
    for i = 2:numel(sample_vec)
        sample_vec(i) = sample_vec(i-1) + fs/sym_rate;
    end
    
    if sample_vec(end) > numel(wave)
        continue;
    end
    
    %%%% Output
    data = [MFout_real(sample_vec); MFout_imag(sample_vec)];
    pack = (sign(reshape(data, 1, n_bits))+1)/2;
    const = data(1,:)+1j*data(2,:);
    eyed = struct('fsfd', fs/sym_rate, 'r', MFout_real(signal_start:sample_vec(end)) + 1j*MFout_imag(signal_start:sample_vec(end)));
    psd = struct('p',[],'f',[]);
    [psd.p,psd.f] = pwelch(wave,[],[],[],fs,'centered','power');
end
stop(rec)
disp('receiver stopped')
end