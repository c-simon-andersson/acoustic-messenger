% clear all
close all
clc

n_bits = 432;
bits = randsrc(n_bits, 1, [0 1]); 
save('bits.mat', 'bits');
transmitter(bits, 2e3);
