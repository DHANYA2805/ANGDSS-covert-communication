function [ber, nErr, nBits] = mc_ber_chain(scheme, N, snr_db, maxBits, minErr, chunkSym)
%MC_BER_CHAIN  Monte Carlo BER of the full QPSK DSSS / GDSS chain (static AWGN).
%   [ber, nErr, nBits] = mc_ber_chain(scheme, N, snr_db, maxBits, minErr, chunkSym)
%   scheme : 'dsss' or 'gdss'
%   Chain  : bits -> QPSK -> PN spreading (long code) -> [Gaussian mask]
%            -> AWGN -> PN despreading -> hard decision -> bit errors
%   SNR is the chip-level Es/N0 of the paper (unit signal power, noise
%   variance 1/SNR). Simulation of each SNR point stops after minErr bit
%   errors or maxBits bits, whichever comes first.
if nargin < 6, chunkSym = max(1, round(2^21 / N)); end
ber   = zeros(size(snr_db));
nErr  = zeros(size(snr_db));
nBits = zeros(size(snr_db));
for k = 1:numel(snr_db)
    nv = 10^(-snr_db(k) / 10);
    while nErr(k) < minErr && nBits(k) < maxBits
        [sym, bits] = qpsk_mod(chunkSym);
        pn    = 2 * randi([0 1], N, chunkSym) - 1;   % fresh code per symbol
        chips = dsss_spread(sym, pn);
        if strcmpi(scheme, 'gdss')
            chips = gdss_mask(chips, 1);
        end
        r     = chips + cgauss(numel(chips), nv);
        z     = despread(r, pn);
        bHat  = qpsk_demod(z);
        nErr(k)  = nErr(k)  + sum(bHat(:) ~= bits(:));
        nBits(k) = nBits(k) + 2 * chunkSym;
    end
    ber(k) = nErr(k) / nBits(k);
end
end
