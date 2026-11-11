function [ber, nErr, nBits] = mc_ber_env(scheme, N, Ns, g, sm2, nv0, minErr, maxBits)
%MC_BER_ENV  Monte Carlo BER through the time-varying noise environment.
%   [ber, nErr, nBits] = mc_ber_env(scheme, N, Ns, g, sm2, nv0, minErr, maxBits)
%   scheme : 'dsss', 'gdss' or 'angdss'
%   N, Ns  : spreading factor, QPSK symbols per block
%   g      : K x 1 noise multipliers (noise_environment), noise var = nv0*g(k)
%   sm2    : K x 1 transmit power per block (ignored for fixed schemes, = 1)
%   The full chain (QPSK, spreading, masking, AWGN, despreading, decision) is
%   simulated block by block; the receiver is identical for all schemes.
%   Blocks are drawn at random from the trace (unbiased average over the
%   noise states) until minErr errors or maxBits bits have been simulated.
K   = numel(g);
B   = max(1, floor(2^20 / (Ns * N)));      % blocks per vectorised chunk
nErr = 0; nBits = 0; nChunk = 0;
while (nErr < minErr || nChunk < 8) && nBits < maxBits   % >= 8 chunks: ~1000 noise states
    nChunk = nChunk + 1;
    idx  = randi(K, 1, B);                 % random noise states
    nSym = B * Ns;
    [sym, bits] = qpsk_mod(nSym);
    pn    = 2 * randi([0 1], N, nSym) - 1;
    chips = dsss_spread(sym, pn);
    perChip = @(v) reshape(repmat(v(idx).', Ns * N, 1), [], 1);
    switch lower(scheme)
        case 'dsss',   tx = chips;
        case 'gdss',   tx = gdss_mask(chips, 1);
        case 'angdss', tx = gdss_mask(chips, sqrt(perChip(sm2)));
        otherwise, error('unknown scheme');
    end
    r    = tx + cgauss(numel(tx), nv0 * perChip(g));
    z    = despread(r, pn);
    bHat = qpsk_demod(z);
    nErr  = nErr + sum(bHat(:) ~= bits(:));
    nBits = nBits + 2 * nSym;
end
ber = nErr / nBits;
end
