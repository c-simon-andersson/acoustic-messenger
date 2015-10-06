function transmitter_v3(packet,fc)
% Complete the function and remove the following lines

sym_rate = 240;
fs = 48e3;
barker = [1 1 1 0 0 0 1 0 0 1 0];
pilot = ones(1, 30);
phaseerror = 0; frequencyerror = fc*0.0001;

% insert barker preamble and pilot for each of the two signals
packet =  [reshape([barker;barker], 1, 2*numel(barker)), reshape([pilot;pilot], 1, 2*numel(pilot)), packet]; 

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
output = real(output).*cos(2*pi*(fc+frequencyerror)*t+phaseerror) + imag(output).*sin(2*pi*(fc+frequencyerror)*t+phaseerror);

% For HIL testing:
% audio_obj = audioplayer(output,fs);
% playblocking(audio_obj);

save('wave.mat', 'output');

end
