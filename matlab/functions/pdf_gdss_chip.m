function p = pdf_gdss_chip(x, sigma_n, sigma_s)
%PDF_GDSS_CHIP  PDF of one noisy GDSS chip r = s + n (base paper Eq. 11).
%   s ~ half-normal(sigma_s), n ~ N(0, sigma_n^2).
%   p(x) = exp(-x^2/(2(ss^2+sn^2))) * (1 + erf(ss*x/(sn*sqrt(2(ss^2+sn^2)))))
%          / sqrt(2*pi*(ss^2+sn^2))
if nargin < 3, sigma_s = 1; end
v = sigma_s^2 + sigma_n^2;
p = exp(-x.^2 / (2*v)) .* (1 + erf(sigma_s * x ./ (sigma_n * sqrt(2*v)))) ./ sqrt(2*pi*v);
end
