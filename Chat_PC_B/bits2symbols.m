function y = bits2symbols(bits, n)
% Converts a vector of bits into symbols with n bits each

symbols = buffer(bits, n)';
s = []; y = 0;
if n == 2    
    s = [(1 + 1i) (1 - 1i) (-1+1i) (-1-1i)] / sqrt(2);
elseif n == 4 %16QAM with Gray labeling
    s = [-3-3i -3-1i -3+3i -3+1i -1-3i -1-1i -1+3i -1+1i 3-3i 3-1i 3+3i 3+1i 1-3i 1-1i 1+3i 1+1i];    
    s = s / sqrt(10);
elseif n == 3 %C8QAM with quasi-grey labeling
    a = sqrt(3);
    s = [-a -a*1i a*1i a -1-1i 1-1i -1+1i 1+1i];
    s = s / sqrt(20/8);
end

symbol_index = bi2de(symbols, 'left-msb')' + 1;
y = s(symbol_index);
end
