function transmitter(packet,fc)
sym_rate = 240;
%sym_rate = 120;
fs = 24e3;
packet = packet';
packet = [packet zeros(1, 3)];
bits_per_sym = 4;
code_length = 1;
barker_norm = [0 sqrt(2) 0 sqrt(10)];
barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1] / sqrt(2);
pilot = ones(1, 30) / sqrt(2);
phase_error = 0;

% Insert barker code for each of the two streams
symbols = [barker+1i*barker (pilot+pilot*1i) bits2symbols(packet, bits_per_sym, code_length)];
symbols_upsampled = upsample(symbols, fs/sym_rate);
%scatter(real(symbols), imag(symbols))

% a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4;
rrc_pulse = rtrcpuls(a,tau,fs,span);

% Construct output waveform
output = conv(symbols_upsampled, rrc_pulse);
t = (1:numel(output))/fs;
output = real(output).*cos(2*pi*fc*t+phase_error) + imag(output).*sin(2*pi*fc*t+phase_error);
output = [zeros(1, 2000) output zeros(1, 2000)];

% Center amplitude around 0 and normalize amplitude to avoid clipping
output = output - mean(output);
output = output / max(abs(output));

audio_obj = audioplayer(output,fs);
disp('transmitting')

output = awgn(output, 10);
save('wave.mat', 'output');
%playblocking(audio_obj);
disp('transmission complete')

end