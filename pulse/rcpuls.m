function y = rcpuls(a, tau, t)

%tau = 1;
t = t + 0.0000001;

tpi = pi/tau; atpi = tpi*a; at = 4*a^2/tau^2;
y = sin(tpi*t) .* cos(atpi*t) ./ (tpi*t .* (1-at*t.^2));

end

