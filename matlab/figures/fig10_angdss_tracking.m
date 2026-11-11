function out = fig10_angdss_tracking(cfg, E)
%FIG10  ANGDSS adaptive loop in a drifting noise floor (sense -> estimate -> adapt).
sa  = cfg.env.demoSNR;                 % average SNR for the demo, e.g. -15 dB
nv0 = 10^(-sa/10);
kk  = 1:min(400, numel(E.g));
gdB  = 10*log10(E.g(kk));
eDB  = 10*log10(E.s2hat(kk) / mean(E.s2hat));
aDB  = 10*log10(E.sm2(kk));
snrG = 10*log10(1 ./ (nv0 * E.g(kk)));
snrA = 10*log10(E.sm2(kk) ./ (nv0 * E.g(kk)));
fh = figure('Visible', cfg.visible, 'Position', [80 80 1000 640]);
subplot(2,1,1); hold on;
plot(kk, gdB, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 4);
plot(kk, eDB, '.', 'Color', [0.20 0.20 0.20], 'MarkerSize', 6);
plot(kk, zeros(size(kk)), '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.8);
plot(kk, aDB, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4);
grid on; box on; xlim([kk(1) kk(end)]);
ylabel('Power relative to mean (dB)');
title('(a) Ambient noise, transmitter estimate and mask power');
legend('True noise floor \sigma_n^2', 'Estimate \sigma_n^2 (median, L samples)', ...
       'GDSS mask \sigma_m^2 (fixed)', 'ANGDSS mask \sigma_m^2 (adaptive)', 'Location', 'southoutside', 'Orientation', 'horizontal');
subplot(2,1,2); hold on;
plot(kk, snrG, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.6);
plot(kk, snrA, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
plot(kk, sa * ones(size(kk)), 'k:', 'LineWidth', 1.2);
grid on; box on; xlim([kk(1) kk(end)]);
xlabel('Block index k'); ylabel('Chip SNR (dB)');
title(sprintf('(b) Instantaneous chip SNR at the receiver/warden, average SNR = %d dB', sa));
legend('GDSS (fixed mask)', 'ANGDSS (noise-matched mask)', 'Design SNR', 'Location', 'southoutside', 'Orientation', 'horizontal');
save_fig(fh, cfg, 'fig10_angdss_tracking');
allG = 10*log10(1 ./ (nv0 * E.g));  allA = 10*log10(E.sm2 ./ (nv0 * E.g));
out.stdSnrGDSS   = std(allG);
out.stdSnrANGDSS = std(allA);
out.estErr_dB    = std(10*log10(E.s2hat ./ E.g));
fprintf('  Fig 10: std of chip SNR: GDSS %.2f dB, ANGDSS %.2f dB; estimator error std %.2f dB\n', ...
    out.stdSnrGDSS, out.stdSnrANGDSS, out.estErr_dB);
end
