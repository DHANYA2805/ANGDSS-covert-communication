function ber = ber_gdss_gauss(N, snr_lin)
%BER_GDSS_GAUSS  Gaussian (CLT) approximation of uncoded QPSK-GDSS BER.
%   Per I/Q branch, sigma_s = 1 and sigma_n^2 = 1/SNR (base paper Eq. 12):
%   z = sum_{i=1}^N (|u_i| + n_i),  E[z] = N*sqrt(2/pi),
%   Var[z] = N*((1 - 2/pi) + 1/SNR),  BER = Q(E[z]/sqrt(Var[z])).
%   At low SNR this tends to DSSS with an SNR loss of 10*log10(pi/2) = 1.96 dB.
mu  = N * sqrt(2/pi);
v   = N .* ((1 - 2/pi) + 1 ./ snr_lin);
ber = 0.5 * erfc(mu ./ sqrt(2 * v));
end
