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
elseif n == 0 %16QAM
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
elseif n == 4 %C16QAM
    % Source: An Optimal Circular 16-QAM Modulation
    % Technique for Wavelet and Fourier Based OFDM
    % - Khaizuran Abdullah, Ahmad Fadzil Ismail, 
    % Wahidah Hashim ∗ and Zahir M. Hussain
    
    samples = samples(1,:) + 1i*samples(2,:);
    
    % Symbol amplitudes
    r1 = sqrt(2);
    r2 = sqrt(3);
    
    P0 = atan(1/r2);
    Ph = pi/3 + P0;
    d = 1;    
    Psi = pi - pi/4 -P0;
    phi = 2*pi - 2*Psi;
    Pp = phi - pi/3;
    b = pi - Pp;
    ds =  sqrt(2*(2*d)^2 - 8*d*cos(b));
    
    r3 = sqrt((1+r2)^2 + 4 - 4*(1+r2)*cos(Ph));
    r4 = sqrt(ds^2 + r1^2 - 2*ds*r1*cos(Pp/2 + Psi));
    
    Vc = d*[r1 1+r2 r3 r4];
    vc = Vc' * [1 1 1 1];
    Avc = vc';
    Avc = reshape(Avc, 1, 16)
    
    % Symbol phases
    t1 = pi/4 * [1 3 5 7 0 2 4 6];
    h0 = asin(2/r3 * sin(Ph));
    h1 = pi - h0;
    h2 = pi + h0;
    h3 = -h0;
    g0 = pi/4 + asin(ds/r4 * sin(Pp/2 + Psi));
    g1 = pi - g0;
    g2 = pi + g0;
    g3 = -g0;
    
    t2 = [h0 h1 h2 h3];
    t3 = [g0 g1 g2 g3];
    
    thetac = [t1 t2 t3];
    
    % Combine ampltiude and phase
    Scir = Avc.*cos(thetac) + 1i*Avc.*sin(thetac);
    s = reshape(Scir, 1, 16);
    s = s / sqrt(mean(abs(s).^2));
    
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

