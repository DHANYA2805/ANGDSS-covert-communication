function out = fig03_gaussianity(cfg)
%FIG03  Distribution and even moments of noise-free GDSS (paper Fig. 7, Table 1).
N      = 256;
nBits  = cfg.tableBits;                % paper: 1e6 bits
nSym   = nBits / 2;
chunk  = 4096;                         % symbols per chunk
orders = 2:2:20;
Sg = []; Sd = [];
xs = [];
done = 0;
while done < nSym
    m = min(chunk, nSym - done);
    [sym, ~] = qpsk_mod(m);
    pn   = 2 * randi([0 1], N, m) - 1;
    d    = dsss_spread(sym, pn);
    g    = gdss_mask(d, 1);
    Sg   = even_moments(Sg, [real(g); imag(g)], orders);
    Sd   = even_moments(Sd, [real(d); imag(d)], orders);
    if isempty(xs), xs = [real(g); imag(g)]; end
    done = done + m;
end
Mg = standardized_moments(Sg);
Md = standardized_moments(Sd);
Mt = arrayfun(@(k) factorial(k) / (2^(k/2) * factorial(k/2)), orders);   % Eq. (4)

% ---- Figure: amplitude histogram vs standard Gaussian ----
z  = (xs - mean(xs)) / std(xs);
ed = linspace(-5, 5, 101);
c  = histc(z, ed); c = c(1:end-1);
ctr = (ed(1:end-1) + ed(2:end)) / 2;
pdfEst = c(:).' / (sum(c) * (ed(2) - ed(1)));
fh = figure('Visible', cfg.visible, 'Position', [100 100 760 480]); hold on;
bar(ctr, pdfEst, 1, 'FaceColor', [0.30 0.60 0.85], 'EdgeColor', 'none');
xx = linspace(-5, 5, 400);
plot(xx, exp(-xx.^2/2)/sqrt(2*pi), 'Color', [0.85 0.33 0.10], 'LineWidth', 2.2);
grid on; box on; xlim([-5 5]);
xlabel('Standardised signal amplitude'); ylabel('Probability density');
title(sprintf('Noise-free GDSS chips vs standard Gaussian (N = %d)', N));
legend('GDSS (simulated)', 'N(0,1) theory', 'Location', 'northeast');
save_fig(fh, cfg, 'fig03_gdss_distribution');

% ---- Table 1 ----
T = [orders(:) Mt(:) Mg(:) Md(:)];
if cfg.saveFigs
    fid = fopen(fullfile(cfg.outdir, 'table1_moments.csv'), 'w');
    fprintf(fid, 'Moment,SGD_theory,GDSS,DSSS\n');
    fprintf(fid, '%d,%.0f,%.0f,%.4f\n', T.');
    fclose(fid);
end
fprintf('  Table 1 (%g bits, N = %d):\n', nBits, N);
fprintf('    %5s %16s %16s %8s\n', 'order', 'Gaussian', 'GDSS', 'DSSS');
fprintf('    %5d %16.0f %16.0f %8.3f\n', T.');
out.table = T;
end
