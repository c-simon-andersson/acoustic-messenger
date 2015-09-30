function y = bits2symbols(bits)
% Converts a vector of bits into symbols with n bits each

y = reshape((bits*2-1), 2, numel(bits)/2);
y = y(1,:) + 1j*y(2,:);
end
