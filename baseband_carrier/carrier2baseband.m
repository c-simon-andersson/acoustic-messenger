function yb = carrier2baseband(pulsetr, fc, t)
%carrier2baseband Shift a pulse train on a carrier down to baseband
%   pulsetr: Pulse train to shift down
%   fc: The carrier frequency
%   t: Vector with time values corresponding to the pule train samples

ybase = pulsetr .* (sqrt(2) .* cos(2*pi*fc*t));

d = fdesign.lowpass('Fp,Fst,Ap,Ast', 0.0001,0.9,0.1,60);
Hd = design(d, 'butter');
yb = filter(Hd, ybase);

end

