function out = fig06_spectral_correlation(cfg)
%FIG06  Cyclostationary (spectral-coherence) detector at -12 dB (paper Fig. 10).
%   Panels: noise only, DSSS with a repeating PN code, GDSS with the same
%   repeating code, and GDSS with a long (per-symbol) code. The x-axis is
%   zoomed to the first 16 PN-period harmonics (alpha = k/256); the same
%   pattern repeats up to alpha = 0.5.
N      = 256;
nSym   = cfg.scfSym;                   % 2048 symbols -> 524,288 chips
Nw     = 4 * N;                        % segment length (multiple of N)
snr    = 10^(-12/10);
L      = N * nSym;
[sym, ~] = qpsk_mod(nSym);
pnS    = 2 * randi([0 1], N, 1) - 1;          % short, repeating code
pnL    = 2 * randi([0 1], N, nSym) - 1;       % long code
dS     = dsss_spread(sym, pnS);
gS     = gdss_mask(dS, 1);
gL     = gdss_mask(dsss_spread(sym, pnL), 1);
sigs   = {zeros(L,1), dS, gS, gL};
names  = {'(a) Gaussian noise', '(b) DSSS, repeating PN, -12 dB', ...
          '(c) GDSS, repeating PN, -12 dB', '(d) GDSS, long PN, -12 dB'};
fh = figure('Visible', cfg.visible, 'Position', [60 60 1250 620]);
out.excess = zeros(1, 4);
for i = 1:4
    x = sqrt(snr) * sigs{i} + cgauss(L, 1);
    [alpha, ees, fl] = spectral_coherence_ees(x, Nw);
    m   = round(alpha * Nw);                                  % cyclic bin index
    on  = (mod(m, Nw / N) == 0) & (m > Nw / N);               % alpha = k/N, k >= 2
    out.excess(i) = mean(ees(on)) / fl - 1;                   % feature strength
    subplot(2, 2, i);
    plot(alpha, ees / fl, 'Color', [0.00 0.45 0.74], 'LineWidth', 1);
    grid on; box on; xlim([0 0.0625]); ylim([0.8 2.4]);
    xlabel('Cyclic frequency \alpha (cycles/chip)'); ylabel('EES / noise floor');
    title(names{i});
end
save_fig(fh, cfg, 'fig06_spectral_correlation');
out.floor = fl;
fprintf('  Fig 6: mean coherence excess at alpha=k/N: noise %.3f, DSSS %.3f, GDSS(short PN) %.3f, GDSS(long PN) %.3f\n', out.excess);
end
