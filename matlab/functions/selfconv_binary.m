function [yN, nConv] = selfconv_binary(y1, N)
%SELFCONV_BINARY  N-fold self-convolution with the paper's fast algorithm.
%   [yN, nConv] = selfconv_binary(y1, N)
%   Builds y2 = y1*y1, y4 = y2*y2, ... (Eq. 17) and then combines only the
%   powers of two present in binary N (Eq. 18). Each '*' is one numerical
%   convolution; nConv returns how many were done (equals Eq. 19).
%   Convolutions are circular (FFT based) on a grid that is wide enough for
%   the sum, see ber_gdss_numeric.m.
y1    = y1(:);
F1    = fft(y1);
bitsN = fliplr(dec2bin(N) == '1');     % LSB first
m     = numel(bitsN) - 1;              % floor(log2 N)
nConv = 0;
P     = cell(m + 1, 1);
P{1}  = F1;                            % spectrum of y_{2^0}
for k = 2:m + 1                        % y_{2^(k-1)} = y_{2^(k-2)} * y_{2^(k-2)}
    P{k}  = P{k-1} .^ 2;
    nConv = nConv + 1;
end
acc = [];
for k = 1:m + 1
    if bitsN(k)
        if isempty(acc)
            acc = P{k};
        else
            acc = acc .* P{k};
            nConv = nConv + 1;
        end
    end
end
yN = real(ifft(acc));
yN(yN < 0) = 0;
yN = yN / sum(yN);
end
