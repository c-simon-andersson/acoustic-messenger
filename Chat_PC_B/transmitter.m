function transmitter(packet,fc)
% Complete the function and remove the following lines

sym_rate = 240;
fs = 24e3;
barker = [1 1 1 0 0 0 1 0 0 1 0];
barker = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
packet = packet';
bits_per_sym = 2;
phase_error = pi/4;

% Insert barker code for each of the two streams
symbols = [(barker+barker*1i) bits2symbols(packet, bits_per_sym)];
symbols_upsampled = upsample(symbols, fs/sym_rate);

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
%playblocking(audio_obj);
disp('transmission complete')

end