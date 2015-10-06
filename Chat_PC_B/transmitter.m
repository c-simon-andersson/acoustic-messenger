function transmitter(packet,fc)
% Complete the function and remove the following lines

sym_rate = 240;
fs = 24e3;
barker = [1 1 1 0 0 0 1 0 0 1 0];
barker = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
packet = packet';
bits_per_sym = 2;
phase_error = pi/4;

% insert barker preamble for each of the two signals
%packet =  [reshape(repmat(barker, bits_per_sym, 1), 1, bits_per_sym*numel(barker)), packet];
%symbols = [bits2symbols(packet(1:length(barker), 2)) bits2symbols(packet(length(barker)+1:end), bits_per_sym)];
symbols = [(barker+barker*1i) bits2symbols(packet, bits_per_sym)]
%symbols = bits2symbols(packet, bits_per_sym);
%scatter(real(symbols), imag(symbols));
symbols_upsampled = upsample(symbols, fs/sym_rate);

% a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4;
rrc_pulse = rtrcpuls(a,tau,fs,span);

output = conv(symbols_upsampled, rrc_pulse);

t = (1:numel(output))/fs;
output = real(output).*cos(2*pi*fc*t+phase_error) + imag(output).*sin(2*pi*fc*t+phase_error);
size(output)

audio_obj = audioplayer(output,fs);
disp('transmitting')
save('wave.mat', 'output');
%playblocking(audio_obj);
disp('transmission complete')

end