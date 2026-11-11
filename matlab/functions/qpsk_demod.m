function bitsHat = qpsk_demod(z)
%QPSK_DEMOD  Hard-decision Gray QPSK demodulation (inverse of qpsk_mod).
%   bitsHat = qpsk_demod(z) returns an nSym x 2 matrix of detected bits.
bitsHat = [real(z(:)) < 0, imag(z(:)) < 0];
end
