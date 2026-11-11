function out = fig04_moment_detector(cfg)
%FIG04  20th-order-moment detector, DSSS vs GDSS (paper Fig. 8).
%   For each SNR the received I/Q samples are standardised (Eq. 5), the 20th
%   moment is estimated and compared with 654,729,075 (Eq. 4). The average
%   absolute deviation in % over nTrials tests is plotted. A noise-only
%   reference is added: it shows the estimator's own spread, because the
%   20th moment of a finite Gaussian sample is itself a noisy statistic
%   (relative std ~ 862/sqrt(n)).
N       = 256;
nBits   = cfg.momentBits;           % paper: 1e6 bits
nTrials = cfg.momentTrials;         % paper: 5
snr_db  = cfg.momentSNR;            % paper: -20 ... 0 dB
M20     = factorial(20) / (2^10 * factorial(10));   % 654,729,075
nSym    = nBits / 2;
chunk   = 4096;
devD = zeros(nTrials, numel(snr_db));
devG = zeros(nTrials, numel(snr_db));
devN = zeros(nTrials, 1);
for t = 1:nTrials
    SD = cell(1, numel(snr_db)); SG = SD; SN = [];
    pn = 2 * randi([0 1], N, 1) - 1;                 % repeating PN code
    done = 0;
    while done < nSym
        m = min(chunk, nSym - done);
        [sym, ~] = qpsk_mod(m);
        d = dsss_spread(sym, pn);
        g = gdss_mask(d, 1);
        L = numel(d);
        for s = 1:numel(snr_db)
            nv = 10^(-snr_db(s)/10);
            rd = d + cgauss(L, nv);
            rg = g + cgauss(L, nv);
            SD{s} = even_moments(SD{s}, [real(rd); imag(rd)], [2 20]);
            SG{s} = even_moments(SG{s}, [real(rg); imag(rg)], [2 20]);
        end
        w  = cgauss(L, 1);
        SN = even_moments(SN, [real(w); imag(w)], [2 20]);
        done = done + m;
    end
    for s = 1:numel(snr_db)
        Md = standardized_moments(SD{s}); Mg = standardized_moments(SG{s});
        devD(t, s) = abs(Md(2) - M20) / M20 * 100;
        devG(t, s) = abs(Mg(2) - M20) / M20 * 100;
    end
    Mn = standardized_moments(SN);
    devN(t) = abs(Mn(2) - M20) / M20 * 100;
    fprintf('    moment test trial %d/%d done\n', t, nTrials);
end
out.snr_db = snr_db;
out.devDSSS  = mean(devD, 1);
out.devGDSS  = mean(devG, 1);
out.devNoise = mean(devN);
% DSSS detection limit: lowest SNR where the deviation exceeds 10 %
i = find(out.devDSSS > 10, 1, 'first');
if isempty(i), out.dsssLimit = NaN; else, out.dsssLimit = snr_db(i); end

% analytic DSSS curve: r = a*c + n, c = +/-1, xi = a^2/sigma^2 (per I/Q branch)
sf = -20:0.25:0; devTh = zeros(size(sf));
dfac = @(k) prod(k:-2:1);                         % double factorial k!!
for s = 1:numel(sf)
    xi = 10^(sf(s)/10); m = 0;
    for j = 0:10
        m = m + nchoosek(20, 2*j) * xi^j * dfac(19 - 2*j);
    end
    devTh(s) = abs(m / (1 + xi)^10 - M20) / M20 * 100;
end
out.devDSSStheory = devTh; out.snrTheory = sf;
i = find(devTh > 10, 1, 'first');
out.dsssLimitTheory = interp1(devTh(i-1:i), sf(i-1:i), 10);
fh = figure('Visible', cfg.visible, 'Position', [100 100 760 480]); hold on;
plot(sf, devTh, '-', 'Color', [0.96 0.70 0.55], 'LineWidth', 1.5);
plot(snr_db, out.devDSSS, '-o', 'LineWidth', 2, 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', 'w');
plot(snr_db, out.devGDSS, '-s', 'LineWidth', 2, 'Color', [0.00 0.45 0.74], 'MarkerFaceColor', 'w');
plot(snr_db, out.devNoise * ones(size(snr_db)), ':', 'LineWidth', 1.8, 'Color', [0.3 0.3 0.3]);
grid on; box on; xlim([min(snr_db) max(snr_db)]); ylim([0 100]);
hline(10, 'r--', '10 % detection threshold');
xlabel('SNR (dB)'); ylabel('Average deviation of 20th moment (%)');
title(sprintf('High-order-moment detector, N = %d, %g bits x %d tests', N, nBits, nTrials));
legend('DSSS (analytic)', 'DSSS (simulated)', 'GDSS (simulated)', 'Noise only (estimator floor)', 'Location', 'northwest');
save_fig(fh, cfg, 'fig04_moment_detector');
fprintf('  Fig 4: DSSS exceeds 10%% from %g dB (analytic crossing %.2f dB); GDSS max deviation %.1f%%; noise-only %.1f%%\n', ...
    out.dsssLimit, out.dsssLimitTheory, max(out.devGDSS), out.devNoise);
end
