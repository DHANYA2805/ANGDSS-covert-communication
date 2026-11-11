function out = fig05_mod_stripping(cfg)
%FIG05  4th-power modulation-stripping detector (paper Fig. 9) + ANGDSS panel.
N     = 256;
nSym  = cfg.modstripSym;              % 10,000 symbols -> 2.56e6 chips
nBins = 72;
[sym, ~] = qpsk_mod(nSym);
pn    = 2 * randi([0 1], N, 1) - 1;
d     = dsss_spread(sym, pn);
g     = gdss_mask(d, 1);
L     = numel(d);
% ANGDSS in a fluctuating noise floor (3 dB std), average SNR = -5 dB
Ns    = 32;  K = ceil(nSym / Ns);
gEnv  = noise_environment(K, cfg.env.sigma_dB, cfg.env.a, cfg.env.clip_dB);
s2h   = estimate_noise_var(gEnv, cfg.env.L, 'median');
sm2   = angdss_mask_power(s2h, 1, 1, cfg.env.peak);
perChip = @(v) local_trim(reshape(repmat(v(:).', Ns * N, 1), [], 1), L);
a     = gdss_mask(d, sqrt(perChip(sm2)));

cases = { cgauss(L, 1),                      '(a) Gaussian noise only';
          d + cgauss(L, 1),                  '(b) DSSS, SNR = 0 dB';
          g + cgauss(L, 1),                  '(c) GDSS, SNR = 0 dB';
          d + cgauss(L, 10^0.5),             '(d) DSSS, SNR = -5 dB';
          g + cgauss(L, 10^0.5),             '(e) GDSS, SNR = -5 dB';
          a + cgauss(L, 10^0.5 * perChip(gEnv)), '(f) ANGDSS, avg SNR = -5 dB, varying noise'};
fh = figure('Visible', cfg.visible, 'Position', [60 60 1200 640]);
out.peakToMean = zeros(1, 6);
for i = 1:6
    [ctr, h] = phase_hist4(cases{i,1}, nBins);
    subplot(2, 3, i);
    plot(ctr, h, '-o', 'Color', [0.1 0.1 0.1], 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'LineWidth', 1);
    grid on; box on; xlim([-pi pi]); ylim([0.75 1.02]);
    set(gca, 'XTick', [-pi -pi/2 0 pi/2 pi], 'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    title(cases{i,2}); xlabel('Phase of -r^4 (rad)'); ylabel('Normalised intensity');
    out.peakToMean(i) = 1 / mean(h);       % 1 for a flat histogram
end
save_fig(fh, cfg, 'fig05_modulation_stripping');
fprintf('  Fig 5: peak/mean of phase histogram  noise %.3f | DSSS0 %.2f GDSS0 %.3f | DSSS-5 %.2f GDSS-5 %.3f | ANGDSS %.3f\n', out.peakToMean);
end

function y = local_trim(x, L)
y = x(1:L);
end
