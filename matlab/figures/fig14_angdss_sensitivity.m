function out = fig14_angdss_sensitivity(cfg, E)
%FIG14  ANGDSS gain vs noise-fluctuation depth, and robustness of the
%   noise estimator (sensing-window length, impulsive interference).
p  = cfg.env;  N = p.N;  t = 1e-5;
sa = -30:0.1:10;
% (a) required average SNR vs fluctuation std
sig = cfg.env.sigmaSweep;
rD = zeros(size(sig)); rG = rD; rA = rD; rI = rD;
for i = 1:numel(sig)
    rng(cfg.seed + 14);          % common random numbers: same drift shape, scaled by sigma
    g   = noise_environment(p.K, sig(i), p.a, p.clip_dB);
    s2h = estimate_noise_var(g, p.L, p.method, p.pImp, p.impGain_dB);
    sm2 = angdss_mask_power(s2h, 1, p.beta, p.peak);
    smg = angdss_mask_power(g,   1, p.beta, p.peak);
    C   = angdss_ber_curves(sa, g, sm2, N, E.tblG);
    Cg  = angdss_ber_curves(sa, g, smg, N, E.tblG);
    rD(i) = required_snr(sa, C.dsss, t);   rG(i) = required_snr(sa, C.gdss, t);
    rA(i) = required_snr(sa, C.angdss, t); rI(i) = required_snr(sa, Cg.angdss, t);
end
% (b) required average SNR vs sensing window L (sigma = p.sigma_dB)
Ls = cfg.env.Lsweep;
rMedImp = zeros(size(Ls)); rMeanImp = rMedImp; rMeanCln = rMedImp;
for i = 1:numel(Ls)
    rng(cfg.seed + 15 + i);
    s2 = estimate_noise_var(E.g, Ls(i), 'median', p.pImp, p.impGain_dB);
    C  = angdss_ber_curves(sa, E.g, angdss_mask_power(s2, 1, p.beta, p.peak), N, E.tblG);
    rMedImp(i) = required_snr(sa, C.angdss, t);
    s2 = estimate_noise_var(E.g, Ls(i), 'mean', p.pImp, p.impGain_dB);
    C  = angdss_ber_curves(sa, E.g, angdss_mask_power(s2, 1, p.beta, p.peak), N, E.tblG);
    rMeanImp(i) = required_snr(sa, C.angdss, t);
    s2 = estimate_noise_var(E.g, Ls(i), 'mean', 0);
    C  = angdss_ber_curves(sa, E.g, angdss_mask_power(s2, 1, p.beta, p.peak), N, E.tblG);
    rMeanCln(i) = required_snr(sa, C.angdss, t);
end
C0 = angdss_ber_curves(sa, E.g, E.sm2g, N, E.tblG);
rGenie = required_snr(sa, C0.angdss, t);  rFix = required_snr(sa, C0.gdss, t);

fh = figure('Visible', cfg.visible, 'Position', [60 60 1250 480]);
subplot(1,2,1); hold on;
plot(sig, rD, '--s', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.8, 'MarkerFaceColor', 'w');
plot(sig, rG, '-o',  'Color', [0.00 0.45 0.74], 'LineWidth', 1.8, 'MarkerFaceColor', 'w');
plot(sig, rA, '-d',  'Color', [0.85 0.33 0.10], 'LineWidth', 2.0, 'MarkerFaceColor', [0.85 0.33 0.10]);
plot(sig, rI, ':',   'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
grid on; box on;
xlabel('Noise-floor fluctuation, std \sigma (dB)'); ylabel('Average SNR for BER = 10^{-5} (dB)');
title(sprintf('(a) Required SNR vs fluctuation depth, N = %d', N));
legend('DSSS', 'GDSS', 'ANGDSS (estimated)', 'ANGDSS (perfect knowledge)', 'Location', 'northwest');
subplot(1,2,2); hold on;
semilogx(Ls, rMeanImp, '-o', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.8, 'MarkerFaceColor', 'w');
semilogx(Ls, rMedImp,  '-d', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.0, 'MarkerFaceColor', [0.85 0.33 0.10]);
semilogx(Ls, rMeanCln, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);
semilogx(Ls, rGenie * ones(size(Ls)), ':', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
semilogx(Ls, rFix * ones(size(Ls)), '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.4);
set(gca, 'XScale', 'log'); grid on; box on; xlim([Ls(1) Ls(end)]);
yl = ylim; ylim([yl(1) yl(2) + 0.45 * diff(yl)]);          % room for the legend
xlabel('Sensing window L (noise samples per block)'); ylabel('Average SNR for BER = 10^{-5} (dB)');
title(sprintf('(b) Estimator robustness, \\sigma = %g dB, impulses p = %.2f, +%d dB', p.sigma_dB, p.pImp, p.impGain_dB));
legend('Mean estimator, with impulses', 'Median estimator, with impulses', 'Mean estimator, clean noise', ...
       'Perfect knowledge', 'Fixed-mask GDSS', 'Location', 'northeast');
save_fig(fh, cfg, 'fig14_angdss_sensitivity');
out.sigma = sig; out.reqDSSS = rD; out.reqGDSS = rG; out.reqANGDSS = rA; out.reqGenie = rI;
out.L = Ls; out.reqMedImp = rMedImp; out.reqMeanImp = rMeanImp; out.reqMeanClean = rMeanCln;
out.reqGenieRef = rGenie; out.reqFixRef = rFix;
fprintf('  Fig 14a: ANGDSS gain over GDSS at sigma = %s dB: %s dB\n', sprintf('%g ', sig), sprintf('%.2f ', rG - rA));
fprintf('  Fig 14b: L = %s -> median %s | mean(imp) %s dB\n', sprintf('%d ', Ls), sprintf('%.2f ', rMedImp), sprintf('%.2f ', rMeanImp));
end
