function transmitter(packet,fc)
%%%% Definitions
sym_rate = 120;
fs = 48e3;
packet = packet';
bits_per_sym = 4;
barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1] / sqrt(2);
pilot = ones(1, 40) / sqrt(2);

%%%% Error simulation
phase_error = 0;
freq_error = 0;
snr = inf;

%%%% Insert barker code for each of the two streams
symbols = [(barker+1i*barker) (pilot+1i*pilot) bits2symbols(packet, bits_per_sym)];
symbols_upsampled = upsample(symbols, fs/sym_rate);

% a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4;
rrc_pulse = rtrcpuls(a,tau,fs,span);

%%%% Construct output waveform
output = conv(symbols_upsampled, rrc_pulse);
t = (1:numel(output))/fs;
output = real(output).*cos(2*pi*(fc+freq_error)*t+phase_error) + imag(output).*sin(2*pi*(fc+freq_error)*t+phase_error);

%%%% Center wave around 0 and normalize amplitude to avoid clipping
output = output - mean(output);
output = output / max(abs(output));

%%%% Add some zeros before and after to prevent popping noises
%output = [zeros(1, 200) output zeros(1, 200)];

audio_obj = audioplayer(output,fs);
disp('transmitting')

output = awgn(output, snr);
save('wave.mat', 'output');
playblocking(audio_obj);
disp('transmission complete')

end