function [iMatched, qMatched] = matchedFilter(iWaveform, qWaveform, pulse)
%matchedFilter Pass the waveform through a matched filter.
%   iWaveform/qWaveform - The waveforms
%   pulse - The pulse to be used for constructing the filter.

% TODO: Is this correct? Doesn't seem to work right...
b = flip(pulse(:));
iMatched = filter(b, 1, iWaveform);
qMatched = filter(b, 1, qWaveform);

end

