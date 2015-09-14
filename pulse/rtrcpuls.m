function y = rtrcpuls(a, tau, t)
%tau = 1;
%Set symbol time
t = t+.0000001;
%Insert offset to prevent NANs
tpi = pi/tau; amtpi = tpi*(1-a); aptpi = tpi*(1 + a);
ac = 4*a/tau; at = 16*a^2/tau^2;
y = (sin(amtpi*t) + (ac*t).*cos(aptpi*t))./(tpi*t.*(1-at*t.^2));
y = y/sqrt(tau);
%Unit energy