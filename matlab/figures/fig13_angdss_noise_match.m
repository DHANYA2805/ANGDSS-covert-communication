function out = fig13_angdss_noise_match(cfg, E)
%FIG13  "Noise-matching" stealth score: KL divergence between the warden's
%   received amplitude distribution and the distribution of the ambient
%   noise alone (same drifting environment), versus average SNR.
%   The warden standardises the whole record (scale-free test), so only the
%   SHAPE of the distribution matters. ANGDSS scales its mask with the noise
%   floor, so signal+noise is a scaled copy of the noise process and the
%   divergence stays at the finite-sample floor; a fixed mask (GDSS) or a
%   fixed constellation (DSSS) changes the shape and is exposed.
N   = cfg.env.N;  Kw = min(cfg.env.Kw, numel(E.g));  nb = cfg.env.nbW;
sa  = -25:2.5:0;
ed  = -7:0.1:7;
perChip = @(v) reshape(repmat(v(1:Kw).', nb, 1), [], 1);
gC  = perChip(E.g);  aC = perChip(E.sm2);
L   = Kw * nb;
pdfOf = @(y) local_pdf([real(y); imag(y)] / sqrt(mean(abs(y).^2)), ed);
kld   = @(p, q) sum(p .* log(p ./ q));
dD = zeros(size(sa)); dG = dD; dA = dD; dF = dD;
for i = 1:numel(sa)
    nv  = 10^(-sa(i)/10);
    sym = qpsk_mod(L / N);
    pn  = 2 * randi([0 1], N, 1) - 1;
    c   = dsss_spread(sym, pn);
    w0  = cgauss(L, nv * gC);                    % noise-only reference record
    q   = pdfOf(w0);
    yD  = c                        + cgauss(L, nv * gC);
    yG  = gdss_mask(c, 1)          + cgauss(L, nv * gC);
    yA  = gdss_mask(c, sqrt(aC))   + cgauss(L, nv * gC);
    yN  = cgauss(L, nv * gC);                    % second noise record (floor)
    dD(i) = kld(pdfOf(yD), q);  dG(i) = kld(pdfOf(yG), q);
    dA(i) = kld(pdfOf(yA), q);  dF(i) = kld(pdfOf(yN), q);
end
fh = figure('Visible', cfg.visible, 'Position', [80 80 820 520]); hold on;
semilogy(sa, dD, '-s', 'Color', [0.47 0.67 0.19], 'LineWidth', 2, 'MarkerFaceColor', 'w');
semilogy(sa, dG, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', 2, 'MarkerFaceColor', 'w');
semilogy(sa, dA, '-d', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.2, 'MarkerFaceColor', [0.85 0.33 0.10]);
semilogy(sa, dF, ':',  'Color', [0.3 0.3 0.3], 'LineWidth', 1.8);
set(gca, 'YScale', 'log'); grid on; box on; xlim([sa(1) sa(end)]);
xlabel('Average SNR (dB)'); ylabel('KL divergence to ambient-noise PDF (nats)');
title(sprintf('Noise-matching stealth score, drifting noise floor (\\sigma = %g dB)', cfg.env.sigma_dB));
legend('DSSS', 'GDSS (fixed mask)', 'ANGDSS (noise-matched mask)', 'Noise vs noise (sampling floor)', 'Location', 'northwest');
save_fig(fh, cfg, 'fig13_angdss_noise_match');
out.sa = sa; out.dsss = dD; out.gdss = dG; out.angdss = dA; out.floor = dF;
j = find(sa == -10, 1);
fprintf('  Fig 13: KLD at -10 dB: DSSS %.2e, GDSS %.2e, ANGDSS %.2e, floor %.2e\n', dD(j), dG(j), dA(j), dF(j));
j = numel(sa);
fprintf('          KLD at   0 dB: DSSS %.2e, GDSS %.2e, ANGDSS %.2e, floor %.2e\n', dD(j), dG(j), dA(j), dF(j));
end

function p = local_pdf(x, ed)
c = histc(x, ed); c = c(1:end-1); c = c(:);
p = (c + 0.5) / sum(c + 0.5);            % add-half smoothing avoids log(0)
end
