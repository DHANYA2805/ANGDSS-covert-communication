function out = fig08_decision_pdf(cfg)
%FIG08  PDF of the GDSS decision variable: numerical vs Monte Carlo (paper Fig. 13).
%   SNR = -10 dB, sigma_s = 1 (per I/Q branch). The Monte Carlo part sums N
%   received chips r_i = |u_i| + n_i, which is exactly what PN despreading
%   produces for one branch of the QPSK-GDSS receiver.
snr_db = -10; snr = 10^(snr_db/10); sn = sqrt(1/snr);
Nlist  = [16 32 64 128 256 512];
nMC    = cfg.pdfTrials;               % 1e5 decision variables per N
edges  = -100:4:800;
ctr    = (edges(1:end-1) + edges(2:end)) / 2;
col    = lines(numel(Nlist));
fh = figure('Visible', cfg.visible, 'Position', [100 100 820 500]); hold on;
hL = zeros(1, numel(Nlist));
out.maxAbsErr = zeros(1, numel(Nlist));
for i = 1:numel(Nlist)
    N = Nlist(i);
    [~, z, pz] = ber_gdss_numeric(N, snr);
    hL(i) = plot(z, pz, '-', 'Color', col(i,:), 'LineWidth', 1.8);
    zs = zeros(nMC, 1); done = 0; B = max(1, floor(2^22 / N));
    while done < nMC
        m = min(B, nMC - done);
        r = abs(randn(N, m)) + sn * randn(N, m);
        zs(done + (1:m)) = sum(r, 1).';
        done = done + m;
    end
    c = histc(zs, edges); c = c(1:end-1);
    pMC = c(:).' / (nMC * (edges(2) - edges(1)));
    plot(ctr, pMC, '.', 'Color', col(i,:), 'MarkerSize', 9);
    pTh = interp1(z, pz, ctr, 'linear', 0);
    out.maxAbsErr(i) = max(abs(pTh - pMC));
end
grid on; box on; xlim([-100 800]);
xlabel('Decision variable z'); ylabel('p_z(z)');
title(sprintf('Decision-variable PDF, SNR = %d dB (lines: numerical, dots: Monte Carlo)', snr_db));
legend(hL, arrayfun(@(n) sprintf('N = %d', n), Nlist, 'UniformOutput', false), 'Location', 'northeast');
save_fig(fh, cfg, 'fig08_decision_pdf');
fprintf('  Fig 8: max |numerical - MC| PDF error per N: %s\n', sprintf('%.1e ', out.maxAbsErr));
end
