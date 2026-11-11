function [alpha, ees, floorLevel] = spectral_coherence_ees(x, Nw)
%SPECTRAL_COHERENCE_EES  Cyclic spectral coherence -> enhanced envelope spectrum.
%   [alpha, ees, floorLevel] = spectral_coherence_ees(x, Nw)
%   Averaged cyclic periodogram estimate of the spectral correlation
%       S_x^alpha(f) = < X_k(f + alpha) X_k*(f) >_k
%   over K = floor(len/Nw) non-overlapping rectangular segments, normalised to
%   the spectral coherence |C(f,alpha)| = |S^alpha(f)| / sqrt(S^0(f+alpha) S^0(f))
%   and averaged over f, which gives the enhanced envelope spectrum (EES)
%   used by the fast spectral-correlation method (Antoni et al., ref. [41] of
%   the base paper). Cyclic frequency alpha is in cycles/chip, step 1/Nw.
%   Choose Nw as a multiple of the spreading factor so that PN-period
%   features (alpha = k/N) fall exactly on the grid.
%   floorLevel = sqrt(pi/(4K)) is the expected EES of white Gaussian noise.
x  = x(:);
K  = floor(numel(x) / Nw);
X  = fft(reshape(x(1:K*Nw), Nw, K), [], 1);      % Nw x K
P  = mean(abs(X).^2, 2);
M  = Nw / 2 - 1;
ees = zeros(M, 1);
for m = 1:M
    Xs  = circshift(X, -m, 1);
    Ps  = circshift(P, -m, 1);
    Sa  = mean(Xs .* conj(X), 2);
    ees(m) = mean(abs(Sa) ./ sqrt(Ps .* P));
end
alpha = (1:M).' / Nw;
floorLevel = sqrt(pi / (4 * K));
end
