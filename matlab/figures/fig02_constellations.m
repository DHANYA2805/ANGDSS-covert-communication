function out = fig02_constellations(cfg)
%FIG02  Noise-free constellations: QPSK-DSSS vs QPSK-GDSS (paper Figs. 3, 5).
N    = 256;
nSym = 40;
[sym, ~] = qpsk_mod(nSym);
pn    = 2 * randi([0 1], N, 1) - 1;
dsss  = dsss_spread(sym, pn);
gdss  = gdss_mask(dsss, 1);
q     = 1 + (real(dsss) < 0) + 2 * (imag(dsss) < 0);      % quadrant 1..4
col   = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19; 0.49 0.18 0.56];
fh = figure('Visible', cfg.visible, 'Position', [100 100 1000 450]);
subplot(1,2,1); hold on;
for c = 1:4
    plot(real(dsss(q==c)), imag(dsss(q==c)), 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', col(c,:), 'MarkerEdgeColor', 'k');
end
axis([-3 3 -3 3]); axis square; grid on; box on;
xlabel('In-phase (I)'); ylabel('Quadrature (Q)');
title('(a) QPSK-DSSS: 4 fixed points');
subplot(1,2,2); hold on;
for c = 1:4
    plot(real(gdss(q==c)), imag(gdss(q==c)), '.', 'MarkerSize', 6, 'Color', col(c,:));
end
axis([-3 3 -3 3]); axis square; grid on; box on;
xlabel('In-phase (I)'); ylabel('Quadrature (Q)');
title('(b) QPSK-GDSS: masked chips, same quadrant');
save_fig(fh, cfg, 'fig02_constellations');
out.quadrantPreserved = all(sign(real(gdss)) == sign(real(dsss)) & sign(imag(gdss)) == sign(imag(dsss)));
fprintf('  Fig 2: every GDSS chip stays in its DSSS quadrant: %d\n', out.quadrantPreserved);
end
