function s2hat = estimate_noise_var(noiseVar, L, method, pImp, impGain_dB)
%ESTIMATE_NOISE_VAR  Transmitter-side ambient noise sensing (ANGDSS step 1-2).
%   s2hat = estimate_noise_var(noiseVar, L, method, pImp, impGain_dB)
%   For every block k the transmitter listens to L complex noise samples of
%   variance noiseVar(k) during a silent sensing slot and estimates it:
%     'mean'   : sigma^2 = mean(|w|^2)            (ML for pure Gaussian noise)
%     'median' : sigma^2 = median(|w|^2) / ln 2   (robust; |w|^2 is exponential)
%     'genie'  : sigma^2 = noiseVar (perfect knowledge, upper bound)
%   pImp, impGain_dB : optional transient impulsive interference in the
%   sensing slot (each sample hit with probability pImp, +impGain_dB power).
%   Vectorised over blocks (processed in chunks to limit memory).
if nargin < 4 || isempty(pImp),       pImp = 0;        end
if nargin < 5 || isempty(impGain_dB), impGain_dB = 20; end
noiseVar = noiseVar(:);
K = numel(noiseVar);
if strcmpi(method, 'genie')
    s2hat = noiseVar;
    return;
end
s2hat = zeros(K, 1);
B = max(1, floor(2^21 / L));                 % blocks per chunk
for k0 = 1:B:K
    idx = k0:min(K, k0 + B - 1);
    v   = noiseVar(idx).';                   % 1 x b
    e   = (randn(L, numel(idx)).^2 + randn(L, numel(idx)).^2) / 2;   % |w|^2 / sigma^2
    if pImp > 0
        hit = rand(L, numel(idx)) < pImp;
        a   = sqrt(10^(impGain_dB/10) / 2);
        ph  = 2*pi*rand(L, numel(idx));                               % circular noise phase
        wI  = sqrt(e) .* cos(ph);
        wQ  = sqrt(e) .* sin(ph);
        wI  = wI + hit .* a .* randn(L, numel(idx));
        wQ  = wQ + hit .* a .* randn(L, numel(idx));
        e   = wI.^2 + wQ.^2;
    end
    e = e .* repmat(v, L, 1);                 % scale to the block's noise power
    switch lower(method)
        case 'mean',   s2hat(idx) = mean(e, 1).';
        case 'median', s2hat(idx) = median(e, 1).' / log(2);
        otherwise, error('unknown estimator %s', method);
    end
end
end
