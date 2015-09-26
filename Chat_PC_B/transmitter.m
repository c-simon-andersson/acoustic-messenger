function transmitter(packet,fc)
% Complete the function and remove the following lines

sym_rate = 240;
fs = 48e3;
barker = [1 1 1 0 0 0 1 0 0 1 0];

% insert barker preamble for each of the two signals
packet =  [reshape([barker;barker], 1, 2*numel(barker)), packet]; 

symbols = bits2symbols(packet);
symbols_upsampled = upsample(symbols, fs/sym_rate);

% a = rolloff, tau = sym time, fs = sampling freq, span = number of sidelobes
a = 0.35; tau = 1/sym_rate; span = 4;
rrc_pulse = rtrcpuls(a,tau,fs,span);

output = conv(symbols_upsampled, rrc_pulse);

figure
subplot(2,1,1)
plot(real(output))
title('transmitter output, real')

subplot(2,1,2)
plot(imag(output))
title('transmitter output, imag')

t = (1:numel(output))/fs;
output = real(output).*cos(2*pi*fc*t) + imag(output).*sin(2*pi*fc*t);

% For HIL testing:
% audio_obj = audioplayer(output,fs);
% playblocking(audio_obj);

save('wave.mat', 'output');

end