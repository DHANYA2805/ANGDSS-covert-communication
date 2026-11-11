function sm2 = angdss_mask_power(s2hat, Pbar, beta, peakFactor)
%ANGDSS_MASK_POWER  Adaptive noise-matched mask variance (proposed ANGDSS).
%   sm2 = angdss_mask_power(s2hat, Pbar, beta, peakFactor)
%   Mask-variance controller (ANGDSS step 3):
%       sigma_m^2[k] = beta * rho * sigma_hat_n^2[k],   rho = Pbar / <sigma_hat_n^2>
%   i.e. the Gaussian mask is kept a fixed ratio rho below (or above) the
%   current noise floor, so the transmitted chips stay matched to the noise.
%   beta       : safety/back-off knob (base value 1; the novelty proposes
%                0.8 ... 1.2 to absorb estimation error)
%   peakFactor : HPA limit, sigma_m^2 <= peakFactor * Pbar (e.g. 4 = +6 dB)
%   Finally the sequence is renormalised so that mean(sigma_m^2) = beta*Pbar:
%   the average transmit power is the SAME as fixed-mask GDSS (fair test).
if nargin < 3 || isempty(beta),       beta = 1;       end
if nargin < 4 || isempty(peakFactor), peakFactor = 4; end
rho = Pbar / mean(s2hat);
sm2 = beta * rho * s2hat(:);
for it = 1:20                               % clip, then restore the budget
    sm2 = min(sm2, peakFactor * Pbar);
    sm2 = sm2 * (beta * Pbar / mean(sm2));
    if max(sm2) <= peakFactor * Pbar * (1 + 1e-9), break; end
end
end
