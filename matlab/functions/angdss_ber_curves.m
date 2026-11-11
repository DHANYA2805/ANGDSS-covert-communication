function C = angdss_ber_curves(sa_db, g, sm2, N, tblG)
%ANGDSS_BER_CURVES  Semi-analytical BER vs AVERAGE SNR in a varying noise floor.
%   C = angdss_ber_curves(sa_db, g, sm2, N, tblG)
%   sa_db : average SNR grid (dB) = mean transmit power / mean noise power
%   g     : K x 1 noise multipliers (mean 1);   sm2 : K x 1 ANGDSS power (mean 1)
%   tblG  : GDSS BER table (snr_db, ber) for this N (from ber_gdss_numeric)
%   For every block the exact (conditional) chip SNR is formed and the block
%   BER is read from the static-AWGN result; the average over all blocks is
%   the BER of the scheme. Returned fields: dsss, gdss, angdss, adsss
%   (adaptive DSSS, i.e. the same power control without Gaussian masking).
%
%   Speed-up: the block chip SNR in dB is (average SNR) + c_k, so the
%   average over K blocks equals a weighted sum over the histogram of c_k
%   (0.01 dB bins). This is exact to within 0.005 dB and ~100x faster.
g = g(:); sm2 = sm2(:);
[cF, wF] = local_hist(-10*log10(g));             % offsets, fixed power
[cA, wA] = local_hist(10*log10(sm2 ./ g));       % offsets, adaptive power
C.sa_db = sa_db;
C.dsss = zeros(size(sa_db)); C.gdss = C.dsss; C.angdss = C.dsss; C.adsss = C.dsss;
for i = 1:numel(sa_db)
    gF = 10.^((sa_db(i) + cF) / 10);
    gA = 10.^((sa_db(i) + cA) / 10);
    C.dsss(i)   = wF.' * ber_dsss_theory(N, gF);
    C.gdss(i)   = wF.' * ber_table_lookup(tblG, gF);
    C.angdss(i) = wA.' * ber_table_lookup(tblG, gA);
    C.adsss(i)  = wA.' * ber_dsss_theory(N, gA);
end
end

function [c, w] = local_hist(x)
q = round(x * 100);                               % 0.01 dB bins
[u, ~, j] = unique(q);
w = accumarray(j(:), 1) / numel(x);
c = u(:) / 100;
end
