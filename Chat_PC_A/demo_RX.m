close all
[pack, psd, const, eyed] = receiver(inf,2e3);
load('../Chat_PC_B/bits.mat');
if ~isempty(pack)
    scatterplot(const')
%     syms = const(1,:)+1j*const(2,:);
%     min_sum_dist = min(abs(syms));
    bit_errors = sum(abs(bits-pack'));
else
    disp('no message received')
end
sprintf('Number of bit errors: %d' ,bit_errors)