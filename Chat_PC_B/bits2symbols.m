function y = bits2symbols(bits, n)
% Converts a vector of bits into symbols with n bits each

%y = reshape((bits*n-1), n, numel(bits)/n);
%y = y(1,:) + 1j*y(2,:);

symbols = buffer(bits, n)';
s = []; y = 0;
if n == 2
    %value_matrix = [0 1;2 3];
    s = [(1 + 1i) (1 - 1i) (-1+1i) (-1-1i)] / sqrt(2);
elseif n == 4
    values = [-3 -1 1 3];
    s = [values+3i values+1i values-1i values-3i] / sqrt(10);
end

symbol_index = bi2de(symbols, 'left-msb')' + 1;
y = s(symbol_index);
end
