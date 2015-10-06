function pack = samples2bits(samples, n)
%samples2bits Takes a vector of samples and converts to bits.

pack = zeros(1, length(samples)*n);
if n == 2
    %s = [(1 + 1i) (1 - 1i) (-1+1i) (-1-1i)] / sqrt(2);
    value_matrix = [2 0;3 1];
    sr = [-1 1] / sqrt(2); si = [1 -1] / sqrt(2);
    s = [sr;si];
    for k = 1:length(samples)
        v = abs(s - repmat(samples(:,k), 1, n));
        [~, vr] = min(v(1,:));
        [~, vi] = min(v(2,:));
        pack((k-1)*2+1:k*2) = de2bi(value_matrix(vi, vr), n, 'left-msb');
    end
elseif n == 4
    value_matrix = [0:3;4:7;8:11;12:15];
    sr = [-3 -1 1 3] / 10;
    si = [3 1 -1 -3] / 10;
    s = [sr;si];
    for k = 1:length(samples)
        v = abs(s - repmat(samples(:,k), 1, n));
        [~, vr] = min(v(1,:));
        [~, vi] = min(v(2,:));
        pack((k-1)*n+1:k*n) = de2bi(value_matrix(vi, vr), n, 'left-msb');
    end    
end

% y = zeros(1, length(samples)*n);
% for k = 1:length(samples)
%     v = abs(s - [samples(:,k) samples(:,k) samples(:,k) samples(:,k)]);
%     [~, vr] = min(-v(2,:));
%     [~, vi] = min(v(1,:));
%     y((k-1)*4+1:k*4) = de2bi(value_matrix(vr, vi), n);
% end

end

