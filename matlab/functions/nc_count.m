function nc = nc_count(N)
%NC_COUNT  Number of numerical convolutions of the fast algorithm (Eq. 19).
%   nc = floor(log2 N) + W(N_bin) - 1, W = number of ones in binary N.
nc = zeros(size(N));
for k = 1:numel(N)
    nc(k) = floor(log2(N(k))) + sum(dec2bin(N(k)) == '1') - 1;
end
end
