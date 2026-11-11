function M = standardized_moments(S)
%STANDARDIZED_MOMENTS  Moments of z = (r - mean)/std (base paper Eq. 5).
%   M = standardized_moments(S), S from even_moments. Order 2 must be in S.
%   The sample mean of these zero-mean signals is O(1/sqrt(n)); with n of
%   order 1e8 its effect on the 20th moment is < 1e-6 relative, so the
%   moments are normalised by the sample variance only.
m   = S.s1 / S.n;
j2  = find(S.orders == 2, 1);
v   = S.sk(j2) / S.n - m^2;
M   = (S.sk / S.n) ./ v.^(S.orders / 2);
end
