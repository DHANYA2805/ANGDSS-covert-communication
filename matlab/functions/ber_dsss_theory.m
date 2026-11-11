function ber = ber_dsss_theory(N, snr_lin)
%BER_DSSS_THEORY  Uncoded QPSK-DSSS BER, base paper Eq. (1).
%   BER = 0.5*erfc( sqrt( N * (Es/N0) / 2 ) ), SNR = Es/N0 at chip level.
ber = 0.5 * erfc(sqrt(N .* snr_lin / 2));
end
