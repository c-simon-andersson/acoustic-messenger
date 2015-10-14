function pack = samples2bits(samples, n)
%samples2bits Takes a vector of samples and converts to bits.

pack = zeros(1, length(samples)*n);
if n == 2 %4QAM
    %s = [(1 + 1i) (1 - 1i) (-1+1i) (-1-1i)] / sqrt(2);
    value_matrix = [2 0;3 1];
    sr = [-1 1] / sqrt(2);
    si = [1 -1] / sqrt(2);
    s = [sr;si];
    
    for k = 1:length(samples)
        v = abs(s - repmat(samples(:,k), 1, n));
        [~, vr] = min(v(1,:));
        [~, vi] = min(v(2,:));
        pack((k-1)*n+1:k*n) = de2bi(value_matrix(vi, vr), n, 'left-msb');
    end
elseif n == 4 %16QAM
    samples = samples(1,:) + 1i*samples(2,:);    
    s = [-3-3i -3-1i -3+3i -3+1i -1-3i -1-1i -1+3i -1+1i 3-3i 3-1i 3+3i 3+1i 1-3i 1-1i 1+3i 1+1i];
    s = s / sqrt(10);
    for k = 1:length(samples)
        [~, v] = min(abs(s - samples(k)));
        pack((k-1)*n+1:k*n) = de2bi(v - 1, n, 'left-msb');
    end
elseif n == 3 %C8QAM with quasi-grey labeling
    samples = samples(1,:) + 1i*samples(2,:);
    a = sqrt(3);
    s = [-a -a*1i a*1i a -1-1i 1-1i -1+1i 1+1i];
    s = s / sqrt(20/8);
    for k = 1:length(samples)
        [~, v] = min(abs(s - samples(k)));
        pack((k-1)*n+1:k*n) = de2bi(v - 1, n, 'left-msb');
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

