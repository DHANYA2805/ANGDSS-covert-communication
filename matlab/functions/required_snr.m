function s_req = required_snr(snr_db, ber, target)
%REQUIRED_SNR  SNR (dB) at which a BER curve crosses the target BER.
%   Interpolates in (dB, log10 BER). Returns NaN if the target is not reached.
ber = ber(:); snr_db = snr_db(:);
ok  = ber > 0;
s_req = NaN;
lb  = log10(ber(ok)); s = snr_db(ok);
i = find(lb(1:end-1) >= log10(target) & lb(2:end) < log10(target), 1);
if ~isempty(i)
    s_req = s(i) + (log10(target) - lb(i)) * (s(i+1) - s(i)) / (lb(i+1) - lb(i));
end
end
