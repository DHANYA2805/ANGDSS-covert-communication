function [sym, bits] = qpsk_mod(nSym)
%QPSK_MOD  Gray-coded unit-power QPSK symbols.
%   [sym, bits] = qpsk_mod(nSym) returns nSym x 1 complex symbols with
%   E|sym|^2 = 1 and the nSym x 2 bit matrix [b_I b_Q].
%   Mapping (Gray): bit 0 -> +1/sqrt(2), bit 1 -> -1/sqrt(2) on each axis
%       00 -> ( 1+1j)/sqrt(2)    01 -> ( 1-1j)/sqrt(2)
%       10 -> (-1+1j)/sqrt(2)    11 -> (-1-1j)/sqrt(2)
bits = randi([0 1], nSym, 2);
sym  = ((1 - 2*bits(:,1)) + 1j*(1 - 2*bits(:,2))) / sqrt(2);
end
