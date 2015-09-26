clear all
close all
clc

n_bits = 416;
bits = randsrc(1, n_bits, [0 1]); 

transmitter(bits, 2e3);

% Simulate channel - add random delays and noise
wave = load('wave.mat');
output = wave.output;
output = [zeros(1,randi(5e3,1,1)), output, zeros(1,randi(5e3,1,1))];
output = awgn(output,20);
save('wave.mat', 'output');

figure
plot(output)
title('channel')

[pack, psd, const, eyed] = receiver_v2(inf,2e3);

scatterplot(const')
syms = const(1,:)+1j*const(2,:);
min_sum_dist = min(abs(syms));
bit_errors = sum(abs(bits-pack));

delete wave.mat