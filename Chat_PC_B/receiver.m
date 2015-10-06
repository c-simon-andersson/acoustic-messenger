function [pack, psd, const, eyed] = receiver(tout,fc)
% Receiver with barker sequence indentification and phase correction

%%%% Definitions
[pack, psd, const, eyed] = deal([]);
barker = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
n_bits = 432;
syms_per_bit = 2;
sym_rate = 240;
fs = 24e3;

barker_threshold = 120;
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
% TODO: Add low-pass filter
LPMF = match_filter;


%rec = audiorecorder(fs, rec_bits, 1);
%record(rec);
pause(0.5)
tic;
%%%% main loop
disp('recieving')
while toc < tout && isempty(pack)

    % Use for HIL
    %wave = getaudiodata(rec, 'int16');
    
    % Use for simulation
    wave = load('wave.mat'); wave = wave.output; wave = wave';   
    
    wave_end = numel(wave);
     
    wave = double(wave(wave_start:end)');
    wave = wave/max(wave);
    t = (1:numel(wave))/fs;
    if( numel(wave) < 1.1*numel(barker_filter) )
        disp('53')         
        continue;       
    end    

    %%%% Shift signal to baseband
    MFout_real = conv(wave.*cos(2*pi*fc*t)*sqrt(2), LPMF);
    MFout_imag = conv(wave.*sin(2*pi*fc*t)*sqrt(2), LPMF);
    
    %%%% Barker Synchronization
    barker_signal_real = fliplr(conv(fliplr(MFout_real), barker_filter, 'same'));    
    barker_signal_imag = fliplr(conv(fliplr(MFout_imag), barker_filter, 'same'));
    barker_signal_sum = barker_signal_real + barker_signal_imag;
    [max_barker, barker_center] = max(abs(barker_signal_sum));    
    
    if max_barker < barker_threshold && ~exist('sample_vec','var')
        wave_start = wave_end;
        disp('78')        
        continue;
    end
    
    signal_start = barker_center + length(barker_upsampled)/2;
    
%     figure; hold on; grid on;
%     plot(barker_signal_sum)
%     plot(barker_center, barker_signal_sum(barker_center), 'ko')
    
    %%%% Phase synchronization
    mes_angle = angle(MFout_real(barker_center) + 1i*MFout_imag(barker_center));
    ref_angle = angle(1+1i);
    phase_error = mes_angle - ref_angle
    MFout_real = conv(wave.*cos(2*pi*fc*t - phase_error)*sqrt(2), LPMF); 
    MFout_imag = conv(wave.*sin(2*pi*fc*t - phase_error)*sqrt(2), LPMF);

    %%%% Sampling
    sample_vec = zeros(1, n_bits/syms_per_bit);

    sample_vec(1) = signal_start;
    for i = 2:numel(sample_vec)
        sample_vec(i) = sample_vec(i-1) + fs/sym_rate;
    end
    
    if sample_vec(end) > numel(wave)
        disp('99')         
        continue;       
    end
    
    %%%% Output    
    data = [MFout_real(sample_vec); MFout_imag(sample_vec)];    
    pack = samples2bits(data, syms_per_bit);     
    const = data(1,:)+1j*data(2,:);
    eyed = struct('fsfd', fs/sym_rate, 'r', MFout_real(signal_start:sample_vec(end)) + 1j*MFout_imag(signal_start:sample_vec(end)));
    psd = struct('p',[],'f',[]);
    %pwelch does not play nice with the gui for some reason
    %[psd.p,psd.f] = pwelch(wave,[],[],[],fs,'centered','power');
end
%stop(rec)
disp('receiver stopped')
end