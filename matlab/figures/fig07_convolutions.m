function out = fig07_convolutions(cfg)
%FIG07  Number of numerical convolutions of the fast algorithm (paper Fig. 12).
Nv = 1:2500;
nc = nc_count(Nv);
ub = 2 * floor(log2(Nv));
fh = figure('Visible', cfg.visible, 'Position', [100 100 760 460]); hold on;
plot(Nv, ub, '.', 'Color', [0.00 0.45 0.74], 'MarkerSize', 6);
plot(Nv, nc, 'x', 'Color', [0.85 0.33 0.10], 'MarkerSize', 3);
grid on; box on;
xlabel('Spreading factor N'); ylabel('Number of convolutions n_c');
title('Fast self-convolution: n_c = \lfloor log_2N \rfloor + W(N_{bin}) - 1');
legend('Upper bound 2\lfloor log_2 N \rfloor', 'n_c (Eq. 19)', 'Location', 'southeast');
save_fig(fh, cfg, 'fig07_num_convolutions');
% check that the implemented algorithm really uses n_c convolutions
[~, ~, ~, ncUsed] = ber_gdss_numeric(50, 1);
out.nc50 = ncUsed;
fprintf('  Fig 7: N = 50 needs %d convolutions (Eq. 19 gives %d, brute force 49)\n', ncUsed, nc_count(50));
end
