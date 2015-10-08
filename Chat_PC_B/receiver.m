function [pack, psd, const, eyed] = receiver(tout,fc)
% Receiver with barker sequence indentification and phase correction

%%%% Definitions
[pack, psd, const, eyed] = deal([]);
%barker = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
pilot = ones(1, 30);
n_bits = 432;
syms_per_bit = 2;
sym_rate = 240;
fs = 24e3;
rec_bits = 24;

barker_threshold = 120;
wave_start = 1;

% a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4; 
rrc_pulse = rtrcpuls(a,tau,fs,span);
match_filter = rrc_pulse;

barker_upsampled = upsample(barker*2-1, fs/sym_rate);
barker_filter = conv(barker_upsampled, rrc_pulse);

%%%% Filter construction
Fpass = (1 + a) / (2*tau);
Fstop = Fpass*1.2;
Apass = 0.01;
Astop = 80;
%LP_filter = firpm(60, [0 cutoff cutoff*1.1 1], [1 1 0 0]); % TODO: Tune filter parameters
%LP_filter = butter(100, cutoff);
%d = fdesign.lowpass('Fp,Fst,Ap,Ast', cutoff,3*cutoff,0.1,60);
filtSpecs = fdesign.lowpass(Fpass, Fstop, Apass, Astop, fs);
Hd = design(filtSpecs, 'ellip');
%LPMF = conv(match_filter, LP_filter); % Combining LP and match-filter
% TODO: Add low-pass filter
LPMF = match_filter;


rec = audiorecorder(fs, rec_bits, 1);
record(rec);
pause(1)
tic;
%%%% main loop
disp('recieving')
while toc < tout && isempty(pack)
    pause(0.2);
    
    % Use for HIL
    %wave = getaudiodata(rec, 'int16');
    wave = getaudiodata(rec, 'double');
    
    % Use for simulation
    %wave = load('wave.mat'); wave = wave.output; wave = wave';   
    
    wave_end = numel(wave);
     
    %wave = double(wave(wave_start:end)');
    wave = wave(wave_start:end)';
    %max(wave)
    wave = wave/max(wave);
    
    t = (1:numel(wave))/fs;
    
%     figure; grid on;
%     plot(wave)
%     title('wave')
    
    % Check that our recording is longer than the barker filter
    if( numel(wave) < 1.1*numel(barker_filter) )
        disp('WARN: Recording shorter than barker sequence.')         
        continue;       
    end
    
    %%%% Shift signal to baseband
    MFout_real = conv(filter(Hd, wave.*cos(2*pi*fc*t)), LPMF);
    MFout_imag = conv(filter(Hd, wave.*sin(2*pi*fc*t)), LPMF);
    %MFout_real = MFout_real / max(MFout_real); MFout_imag = MFout_imag / max(MFout_imag);
    
%     %%%% Automatic gain control
%     gain = 1 / sqrt(mean(MFout_real.^2 + MFout_imag.^2))
%     MFout_real = MFout_real*gain; MFout_imag = MFout_imag*gain;
    
    %%%% Barker Synchronization
    barker_signal_real = fliplr(conv(fliplr(MFout_real), barker_filter, 'same'));
    barker_signal_imag = fliplr(conv(fliplr(MFout_imag), barker_filter, 'same'));
    barker_signal_sum = barker_signal_real + barker_signal_imag;
    [max_barker, barker_center] = max(abs(barker_signal_sum))    
    
%      figure; grid on;
%      plot(barker_signal_sum)
    
    if max_barker < barker_threshold && ~exist('sample_vec','var')
        wave_start = wave_end;
        disp('INFO: Ran loop without finding barker sequence.')
        continue;
    end
    
    barker_start = barker_center - length(barker_upsampled)/2;
    barker_end = barker_center + length(barker_upsampled)/2;
    signal_start = barker_end + length(pilot)*fs/sym_rate;
    
%     figure; hold on; grid on;
%     plot(barker_signal_sum)
%     plot(barker_center, barker_signal_sum(barker_center), 'ko')

    %%%% Prepare sampling
    sample_vec = zeros(1, n_bits/syms_per_bit);

    sample_vec(1) = signal_start;
    for i = 2:numel(sample_vec)
        sample_vec(i) = sample_vec(i-1) + fs/sym_rate;
    end
    
    if sample_vec(end) > numel(wave)
        disp('INFO: Packet continues in next recording.')         
        continue;       
    end
    
    %%%% Phase synchronization
    %mes_angle = angle(MFout_real(barker_start:fs/sym_rate:signal_start-fs/sym_rate) + 1i*MFout_imag(barker_start:fs/sym_rate:signal_start-fs/sym_rate));
    %mes_angle = angle(MFout_real(barker_end:fs/sym_rate:signal_start-fs/sym_rate) + 1i*MFout_imag(barker_end:fs/sym_rate:signal_start-fs/sym_rate));
    mes_angle = angle(MFout_real(barker_end:signal_start-fs/sym_rate) + 1i*MFout_imag(barker_end:signal_start-fs/sym_rate));
    %ref_angle = angle([barker pilot] + 1i*[barker pilot]);
    %ref_angle = angle([pilot] + 1i*[pilot]);
    ref_angle = angle(1+1i);
    phase_error = mean(mes_angle - ref_angle)
%     ref_angle = angle(1+1i);
%     phase_error = mes_angle - ref_angle
    MFout_real = conv(filter(Hd, wave.*cos(2*pi*fc*t - phase_error)), LPMF); 
    MFout_imag = conv(filter(Hd, wave.*sin(2*pi*fc*t - phase_error)), LPMF);

    %%%% Automatic gain control
    gain = 1 / sqrt(mean(MFout_real.^2 + MFout_imag.^2))
    MFout_real = MFout_real*gain; MFout_imag = MFout_imag*gain;   
    
    %%%% Output    
    data = [MFout_real(sample_vec); MFout_imag(sample_vec)];
    pack = samples2bits(data, syms_per_bit);
    const = data(1,:)+1j*data(2,:);
    eyed = struct('fsfd', fs/sym_rate, 'r', MFout_real(signal_start:sample_vec(end)) + 1j*MFout_imag(signal_start:sample_vec(end)));
    psd = struct('p',[],'f',[]);
    %pwelch does not play nice with the gui for some reason
    %[psd.p,psd.f] = pwelch(wave,[],[],[],fs,'centered','power');
    
    figure; grid on; hold on;
    plot(MFout_real); plot(MFout_imag);
end
stop(rec)
disp('receiver stopped')
end