function [pack, psd, const, eyed] = receiver(tout,fc)
%Receiver with barker sequence indentification, phase  and frequency correction

%%%% Definitions
[pack, psd, const, eyed] = deal([]);
barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1] / sqrt(2);
pilot = ones(1, 40) / sqrt(2);
n_bits = 432;
syms_per_bit = 3;
sym_rate = 120;
fs = 48e3;
rec_bits = 24;

barker_threshold = 120;
wave_start = 1;

% Activate debug plots
debug_plots = 1;

% a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4; 
rrc_pulse = rtrcpuls(a,tau,fs,span);

barker_upsampled = upsample(barker*2-1, fs/sym_rate);
barker_filter = conv(barker_upsampled, rrc_pulse);

%%%% Filter construction
match_filter = rrc_pulse;
LPMF = match_filter;

rec = audiorecorder(fs, rec_bits, 1);
record(rec);
pause(1)
tic;

%%%% main loop
disp('Recieving...')
while toc < tout && isempty(pack)
    pause(0.2);
    
    % Use for HIL    
    %wave = getaudiodata(rec, 'double');
    
    % Use for simulation
    wave = load('wave.mat'); wave = wave.output; wave = wave';   
    
    wave_end = numel(wave);      
    wave = wave(wave_start:end)';   
    barker_threshold = 100*max(wave) + 1;
    
    t = (1:numel(wave))/fs;
    
    % Check that our recording is longer than the barker filter
    if( numel(wave) < 1.1*numel(barker_filter) )
        disp('WARN: Recording shorter than barker sequence.')         
        continue;       
    end
    
    %%%% Shift signal to baseband    
    MFout_real = conv(wave.*cos(2*pi*fc*t), LPMF);
    MFout_imag = conv(wave.*sin(2*pi*fc*t), LPMF);    
    
    %%%% Barker Synchronization
    barker_signal_real = fliplr(conv(fliplr(MFout_real), barker_filter, 'same'));
    barker_signal_imag = fliplr(conv(fliplr(MFout_imag), barker_filter, 'same'));
    barker_signal_sum = sqrt(barker_signal_real.^2 + barker_signal_imag.^2);
    [max_barker, barker_center] = max(abs(barker_signal_sum))

    if max_barker < barker_threshold && ~exist('sample_vec','var')
        wave_start = wave_end;
        disp('INFO: Ran loop without finding barker sequence.')
        continue;
    end
    
    %%%% Use barker to define starting points
    barker_start = barker_center - length(barker_upsampled)/2;
    barker_end = barker_center + length(barker_upsampled)/2;
    signal_start = barker_end + length(pilot)*fs/sym_rate;
    signal_end = signal_start + n_bits/syms_per_bit*fs/sym_rate;

    %%%% Prepare sampling
    sample_vec = zeros(1, n_bits/syms_per_bit);

    sample_vec(1) = signal_start;
    for i = 2:numel(sample_vec)
        sample_vec(i) = sample_vec(i-1) + fs/sym_rate;
    end
    
    sample_vec(end)
    numel(wave)
    if sample_vec(end) > numel(wave)
        disp('INFO: Packet continues in next recording.')
        continue;
    end
    
    %%%% Frequency synchronization
    use_frequency_synch = 0;
    if use_frequency_synch
        mes_freq = angle(MFout_real(barker_end:fs/sym_rate:signal_start-fs/sym_rate) + 1i*MFout_imag(barker_end:fs/sym_rate:signal_start-fs/sym_rate));
        diff_freq = mes_freq(1:end-1) - mes_freq(2:end);
        freq_error = (mean(diff_freq) / (2*pi) * sym_rate)
        fc = fc + freq_error;
        MFout_real = conv(wave.*cos(2*pi*fc*t), LPMF);
        MFout_imag = conv(wave.*sin(2*pi*fc*t), LPMF);
    end
    
    %%%% Phase synchronization    
    mes_angle = angle(MFout_real(barker_end:fs/sym_rate:signal_start-fs/sym_rate) + 1i*MFout_imag(barker_end:fs/sym_rate:signal_start-fs/sym_rate));
    ref_angle = angle(1+1i);
    diff_angle = mes_angle - ref_angle;
    phase_error = mean(diff_angle)
    MFout_real = conv(wave.*cos(2*pi*fc*t - phase_error), LPMF); 
    MFout_imag = conv(wave.*sin(2*pi*fc*t - phase_error), LPMF);

    %%%% Automatic gain control
    gain = 1 / sqrt(mean(MFout_real(signal_start:signal_end).^2 + MFout_imag(signal_start:signal_end).^2))
    MFout_real = MFout_real*gain; MFout_imag = MFout_imag*gain;   
    
    %%%% Output    
    data = [MFout_real(sample_vec); MFout_imag(sample_vec)];
    pack = samples2bits(data, syms_per_bit);
    const = data(1,:)+1j*data(2,:);
    eyed = struct('fsfd', fs/sym_rate, 'r', MFout_real(signal_start:sample_vec(end)) + 1j*MFout_imag(signal_start:sample_vec(end)));
    psd = struct('p',[],'f',[]);
    [psd.p, psd.f] = pwelch(wave, [], [], (fc-400):(fc+400), fs, 'twosided');
    psd.p = 10*log10(psd.p./max(psd.p));
    psd.f = psd.f - fc;
    
    if debug_plots
        figure; grid on; hold on;
        plot(MFout_real); plot(MFout_imag);
        plot(signal_start, MFout_real(signal_start), 'ko')
        plot(signal_end, MFout_real(signal_end), 'ko')        
        
        figure; grid on; hold on;
        plot(barker_signal_sum);
        plot(barker_center, barker_signal_sum(barker_center), 'ko')
    end
end
stop(rec)
disp('Receiver stopped!')
end