function [ber, z, pz, nConv] = ber_gdss_numeric(N, snr_lin)
%BER_GDSS_NUMERIC  Semi-analytical uncoded QPSK-GDSS BER (base paper Sec. 5.1).
%   [ber, z, pz, nConv] = ber_gdss_numeric(N, snr_lin)
%   1) single noisy chip PDF p_r(x), Eq. (11)/(13) with sigma_s = 1,
%      sigma_n^2 = 1/SNR (per I/Q branch, Eq. 12)
%   2) normalised y1 = p_r * dx (Eq. 15)
%   3) N-fold self-convolution with the fast binary algorithm (Eq. 16-19)
%   4) BER = sum of p_z over z < 0 (Eq. 7, 21)
%
%   Improvement over the paper's implementation: the convolutions are done on
%   a circular grid whose width (+/-14 std of z) and step (1/16 of the chip
%   PDF scale) are chosen automatically from N and SNR, so no trial-and-error
%   x_max is needed and the whole Figure 14 takes seconds instead of hours.
%   z and pz are returned for the last SNR value (used for Figure 13).
ber = zeros(size(snr_lin));
for k = 1:numel(snr_lin)
    sn  = sqrt(1 / snr_lin(k));
    mu1 = sqrt(2/pi);
    v1  = (1 - 2/pi) + sn^2;
    muZ = N * mu1;
    sdZ = sqrt(N * v1);
    sr  = sqrt(1 + sn^2);
    dx0 = min(sr, sn * sr) / 16;            % resolve the single-chip PDF
    W   = 28 * sdZ;                          % circular period (+/-14 sigma)
    Lg  = 2^nextpow2(ceil(W / dx0));
    Lg  = min(max(Lg, 2048), 2^21);
    dx  = W / Lg;
    x1  = mu1 + dx * ((0:Lg-1).' - Lg/2);    % chip grid centred on E[r]
    y1  = pdf_gdss_chip(x1, sn) * dx;        % Eq. (15)
    y1  = y1 / sum(y1);
    [yN, nConv] = selfconv_binary(y1, N);     % Eq. (16)-(19)
    z0  = N * x1(1) + dx * (0:Lg-1).';       % values before unwrapping
    zz  = z0 - W * floor((z0 - (muZ - W/2)) / W);   % unwrap around E[z]
    ber(k) = sum(yN(zz < 0));                 % Eq. (21)
end
[z, idx] = sort(zz);
pz = yN(idx) / dx;                            % Eq. (20)
end
