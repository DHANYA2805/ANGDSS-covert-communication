function x = gdss_mask(chips, sigma_m)
%GDSS_MASK  Gaussian masking of spread QPSK chips (base paper, Section 3).
%   x = gdss_mask(chips, sigma_m)
%   I = Re(S).*|U|,  Q = Im(S).*|V|,  U,V ~ N(0, sigma_m^2) independent.
%   sigma_m : scalar (fixed mask = original GDSS) or a vector with one value
%             per chip (block-adaptive mask = proposed ANGDSS).
%   Only the magnitude of each chip is changed, so the sign (quadrant) still
%   carries the data and the receiver despreads exactly as in DSSS.
%   With unit-power QPSK chips and sigma_m = 1, E|x|^2 = 1 (same as DSSS).
U = sigma_m .* abs(randn(size(chips)));
V = sigma_m .* abs(randn(size(chips)));
x = real(chips) .* U + 1j * imag(chips) .* V;
end
