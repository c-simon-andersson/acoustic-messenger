function transmitter(packet,fc)
sym_rate = 240;
sym_rate = 120;
fs = 12e3;
packet = packet';
bits_per_sym = 4;
barker_norm = [0 sqrt(2) 0 3*sqrt(10)];
barker = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];% / barker_norm(bits_per_sym);
phase_error = pi;

% Insert barker code for each of the two streams
symbols = [(barker+barker*1i) bits2symbols(packet, bits_per_sym)];
symbols_upsampled = upsample(symbols, fs/sym_rate);
%scatter(real(symbols), imag(symbols))

% a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4;
rrc_pulse = rtrcpuls(a,tau,fs,span);

% Construct output waveform
output = conv(symbols_upsampled, rrc_pulse);
t = (1:numel(output))/fs;
output = real(output).*cos(2*pi*fc*t+phase_error) + imag(output).*sin(2*pi*fc*t+phase_error);

audio_obj = audioplayer(output,fs);
disp('transmitting')
save('wave.mat', 'output');
playblocking(audio_obj);
disp('transmission complete')

end