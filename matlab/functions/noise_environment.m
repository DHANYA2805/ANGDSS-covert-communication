function g = noise_environment(K, sigma_dB, a, clip_dB)
%NOISE_ENVIRONMENT  Time-varying ambient noise floor (ANGDSS test scenario).
%   g = noise_environment(K, sigma_dB, a, clip_dB) returns K x 1 positive
%   multipliers of the average noise power, one per transmission block.
%   The noise level in dB follows a stationary AR(1) (Gauss-Markov) process
%       x_k = a x_{k-1} + sqrt(1-a^2) * sigma_dB * e_k,  e_k ~ N(0,1)
%   clipped to +/- clip_dB, i.e. a slowly drifting, log-normal noise floor
%   (interference, temperature, other users). g is normalised to mean(g) = 1
%   so that the AVERAGE SNR is exactly the value set by the caller.
if nargin < 3 || isempty(a),       a = 0.95; end
if nargin < 4 || isempty(clip_dB), clip_dB = 10; end
x    = zeros(K, 1);
x(1) = sigma_dB * randn;
for k = 2:K
    x(k) = a * x(k-1) + sqrt(1 - a^2) * sigma_dB * randn;
end
x = max(min(x, clip_dB), -clip_dB);
g = 10.^(x / 10);
g = g / mean(g);
end
