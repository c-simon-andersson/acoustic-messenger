function y = rcpuls(a, tau, t)
t = t + 0.0000001;

tpi = pi/tau; atpi = tpi*a; at = 4*a^2/tau^2;
y = sin(tpi*t) .* cos(atpi*t) ./ (tpi*t .* (1-at*t.^2));
norm = sqrt(sum(y.^2));
y = y/norm; % Unit energy

end

