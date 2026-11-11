function out = fig11_angdss_ber(cfg, E)
%FIG11  BER vs average SNR in a drifting noise floor: DSSS vs GDSS vs ANGDSS.
%   Same average transmit power for every scheme. Lines: semi-analytical
%   average of the exact block BER; markers: Monte Carlo of the full chain.
N  = cfg.env.N;
sa = -30:0.1:5;
C  = angdss_ber_curves(sa, E.g, E.sm2,  N, E.tblG);
Cg = angdss_ber_curves(sa, E.g, E.sm2g, N, E.tblG);
st.gdss = ber_table_lookup(E.tblG, 10.^(sa/10));            % static AWGN reference
fh = figure('Visible', cfg.visible, 'Position', [80 80 860 560]); hold on;
h(1) = semilogy(sa, C.dsss,   '--', 'Color', [0.47 0.67 0.19], 'LineWidth', 2);
h(2) = semilogy(sa, C.gdss,   '-',  'Color', [0.00 0.45 0.74], 'LineWidth', 2);
h(3) = semilogy(sa, C.angdss, '-',  'Color', [0.85 0.33 0.10], 'LineWidth', 2.4);
h(4) = semilogy(sa, Cg.angdss, ':', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
h(5) = semilogy(sa, st.gdss,  '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
% Monte Carlo check
sMC = cfg.env.mcSNR;
nv  = 10.^(-sMC/10);
mcD = zeros(size(sMC)); mcG = mcD; mcA = mcD;
for i = 1:numel(sMC)
    mcD(i) = mc_ber_env('dsss',   N, cfg.env.Ns, E.g, E.sm2, nv(i), cfg.mcMinErr, cfg.env.mcMaxBits);
    mcG(i) = mc_ber_env('gdss',   N, cfg.env.Ns, E.g, E.sm2, nv(i), cfg.mcMinErr, cfg.env.mcMaxBits);
    mcA(i) = mc_ber_env('angdss', N, cfg.env.Ns, E.g, E.sm2, nv(i), cfg.mcMinErr, cfg.env.mcMaxBits);
end
semilogy(sMC, max(mcD,1e-9), 's', 'Color', [0.47 0.67 0.19], 'MarkerSize', 7, 'LineWidth', 1.2);
semilogy(sMC, max(mcG,1e-9), 'o', 'Color', [0.00 0.45 0.74], 'MarkerSize', 7, 'LineWidth', 1.2);
semilogy(sMC, max(mcA,1e-9), 'd', 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 6);
set(gca, 'YScale', 'log'); grid on; box on; ylim([1e-6 0.5]); xlim([-25 5]);
xlabel('Average SNR (dB)'); ylabel('Bit error rate');
title(sprintf('Drifting noise floor (\\sigma = %g dB), N = %d, equal average power', cfg.env.sigma_dB, N));
legend(h, 'DSSS (fixed power)', 'GDSS (fixed mask)', 'ANGDSS (estimated noise)', ...
    'ANGDSS (perfect noise knowledge)', 'GDSS, constant noise (reference)', 'Location', 'southwest');
save_fig(fh, cfg, 'fig11_angdss_ber');
t = 1e-5;
out.req.dsss   = required_snr(sa, C.dsss, t);
out.req.gdss   = required_snr(sa, C.gdss, t);
out.req.angdss = required_snr(sa, C.angdss, t);
out.req.genie  = required_snr(sa, Cg.angdss, t);
out.req.adsss  = required_snr(sa, C.adsss, t);
out.req.static = required_snr(sa, st.gdss, t);
out.mc = struct('snr', sMC, 'dsss', mcD, 'gdss', mcG, 'angdss', mcA);
out.curves = C; out.curvesGenie = Cg;
mcTh = interp1(sa, C.angdss, sMC);
fprintf('  Fig 11: SNR for BER 1e-5: DSSS %.2f | GDSS %.2f | ANGDSS %.2f (genie %.2f) | static GDSS %.2f | adaptive DSSS %.2f dB\n', ...
    out.req.dsss, out.req.gdss, out.req.angdss, out.req.genie, out.req.static, out.req.adsss);
fprintf('          ANGDSS gain over GDSS %.2f dB, over DSSS %.2f dB\n', out.req.gdss - out.req.angdss, out.req.dsss - out.req.angdss);
fprintf('          MC vs theory (ANGDSS): %s\n', sprintf('[%.2e/%.2e] ', [mcA; mcTh]));
end
