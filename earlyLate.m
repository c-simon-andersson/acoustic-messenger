function a = earlyLate(waveform, n, ns, epsilon)
%earlyLate Early late gate synchronization
%   waveform - Waveform to synchronise
%   n - Samples between points
%   ns - Step size
%   Returns the first point on the waveform for which the value difference
%   between the two points is less than epsilon*max(waveform).

epsilon = epsilon*max(waveform);
d = 0;
a = 1; b = a + n;

% Raise waveform above 0 to avoid positive and negative samples canceling
minimum = min(waveform);
if minimum < 0
    waveform = waveform + abs(minimum);
end

% Skip all zero values at the beginning of the waveform
while abs(d) < epsilon
    a = a + ns;
    b = b + ns;
    
    if b > length(waveform)
        a = 0; b = 0;
        break;
    end
    
    d = waveform(a) - waveform(b);
end

% Attempt to find an extreme point
while abs(d) > epsilon
   a = a + ns;
   b = b + ns;
   
   if b > length(waveform)
      a = 0; b = 0;
      break;
   end
   
   d = waveform(a) - waveform(b);
end

a = floor((a + b) / 2);


end

