clear; close; clc;
%bits = [0 0 1 1 0 0 0 1 1 1 0 0];
bits = [1 1 0 0 1 1];
symbols = [0 1];
SNR = 1;
fc = 1000;

%TODO: Implement symbol coding
data = bits; 

alpha = 0.3;
Ts = 1/240;
n = 20000;
width = length(data)*2;
pulse = 'rtrcpuls';
%symbolpulse = feval(pulse, alpha, -width:1/n:width);
symbolpulse = pulsetr(pulse, alpha, Ts, n, width, symbols(1));
%symbolpulses = [symbols(1)*symbolpulse;symbols(2)*symbolpulse];
symbolpulses = [pulsetr(pulse, alpha, Ts, n, 1, symbols(1));
                pulsetr(pulse, alpha, Ts, n, 1, symbols(2))];

y = pulsetr(pulse, alpha, Ts, n, width, data);
%t = [-width/2:1/n:(width/2+length(data)-1)];
t = -Ts*width/2:Ts/n:Ts*width/2+(length(data)-1)*Ts;

%Add AWGN noise
%y = awgn(y, SNR);

% Plot waveform
figure
grid on
hold on
plot(t, y)
plot(0:Ts:(length(data)-1)*Ts, data, 'ko')
xlabel('Time (s)');
ylabel('Amplitude')

% To/from carrier
%d = fdesign.lowpass('Fp,Fst,Ap,Ast', 0.0001,0.0009,0.1,60);
d = fdesign.lowpass('Fp,Fst,Ap,Ast', 1/n,2/n,0.1,60);
Hd = design(d, 'butter');
yc = baseband2coscarrier(y, fc, t);
yb = coscarrier2baseband(yc, fc+fc*0.01, t, Hd);

% Plot waveform at carrier and baseband
figure
grid on
hold on
plot(t, yc)
plot(t, y, 'color', 'k')
plot(t, yb, 'LineWidth', 2)
legend('Signal on carrier', 'Original waveform + awgn with SNR=10', 'After shift to baseband and LPF');
xlabel('Time (s)');
ylabel('Amplitude')

% Interpret data using various receivers
%samplingRxData = samplingrx(yb, symbols, n, width)
%correlatingRxData = correlationrx(y, symbols, symbolpulses, n, width)
%mfrWaveform = matchedfilterrx(y, symbolpulse, n, width);

% FFT of carrier and baseband waveform
NFFT = 2^nextpow2(length(yc)); % Next power of 2 from length of yc
YC = fft(yc,NFFT)/length(yc);
YB = fft(yb, NFFT)/length(yb);
f = n/Ts/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum
figure
hold on
grid on

%plot(f,2*abs(YC(1:NFFT/2+1)))
%plot(f,2*abs(YB(1:NFFT/2+1)))

[~,b] = max(YC);
plot(f(1:2*b),2*abs(YC(1:2*b)))
plot(f(1:2*b),2*abs(YB(1:2*b)))
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

% Plot frequency mask
maxYC = 2*abs(max(YC));
plot([fc-200 fc-150 fc-150 fc-120 fc-120 fc+120 fc+120 fc+150 fc+150 fc+200], [db2mag(-25)*maxYC db2mag(-25)*maxYC db2mag(-10)*maxYC db2mag(-10)*maxYC, maxYC, maxYC, db2mag(-10)*maxYC db2mag(-10)*maxYC db2mag(-25)*maxYC db2mag(-25)*maxYC], 'color', 'g')

%%
oneToZero = sum(data) - sum(data .* samplingRxData);
zeroToOne = sum(samplingRxData) - sum(data .* samplingRxData);

disp('Total flips: ')
disp(oneToZero + zeroToOne)

disp('1 -> 0: ')
disp(oneToZero)

disp('0 -> 1: ')
disp(zeroToOne)


