function out = fig09_ber_uncoded(cfg)
%FIG09  Uncoded BER: QPSK-GDSS (semi-analytical) vs QPSK-DSSS (paper Fig. 14).
%   Solid: GDSS (numerical convolution, Sec. 5.1); dashed: DSSS (Eq. 1);
%   markers: Monte Carlo simulation of the full transmitter/receiver chain.
snr_db = -35:0.25:10;
snr    = 10.^(snr_db/10);
Nlist  = [16 32 64 128 256 512 1024 2048];
col    = [0.80 0.10 0.10; 0.90 0.40 0.05; 0.85 0.65 0.00; 0.40 0.70 0.20; ...
          0.10 0.60 0.55; 0.10 0.45 0.80; 0.35 0.25 0.75; 0.55 0.15 0.55];   % one colour per N
fh = figure('Visible', cfg.visible, 'Position', [80 80 900 560]); hold on;
out.N = Nlist; out.snrG = zeros(size(Nlist)); out.snrD = out.snrG;
hG = zeros(size(Nlist));
for i = 1:numel(Nlist)
    N  = Nlist(i);
    bG = ber_gdss_numeric(N, snr);
    bD = ber_dsss_theory(N, snr);
    hG(i) = semilogy(snr_db, max(bG, 1e-300), '-',  'Color', col(i,:), 'LineWidth', 1.8);
    semilogy(snr_db, max(bD, 1e-300), '--', 'Color', col(i,:), 'LineWidth', 1.2);
    out.snrG(i) = required_snr(snr_db, bG, 1e-5);
    out.snrD(i) = required_snr(snr_db, bD, 1e-5);
    out.curves{i} = bG;
end
out.penalty = out.snrG - out.snrD;
% Monte Carlo points (full chain)
out.mc = struct('N', {}, 'snr_db', {}, 'berG', {}, 'berD', {});
for N = cfg.mcN
    i  = find(Nlist == N);
    sAll = snr_db(1:8:end);                              % every 2 dB
    bGp  = interp1(snr_db, out.curves{i}, sAll);
    bDp  = ber_dsss_theory(N, 10.^(sAll/10));
    sG   = sAll(bGp > cfg.mcMinBER & bGp < 0.4);         % points worth simulating
    sD   = sAll(bDp > cfg.mcMinBER & bDp < 0.4);
    berG = mc_ber_chain('gdss', N, sG, cfg.mcMaxBits, cfg.mcMinErr);
    berD = mc_ber_chain('dsss', N, sD, cfg.mcMaxBits, cfg.mcMinErr);
    semilogy(sG, max(berG, 1e-9), 'o', 'Color', col(i,:), 'MarkerFaceColor', col(i,:), 'MarkerSize', 5);
    semilogy(sD, max(berD, 1e-9), 's', 'Color', col(i,:), 'MarkerSize', 6, 'LineWidth', 1.2);
    out.mc(end+1) = struct('N', N, 'snr_db', sG, 'berG', berG, 'berD', berD);
    out.mcD{numel(out.mc)} = sD;
end
set(gca, 'YScale', 'log'); grid on; box on; ylim([1e-7 1]); xlim([-35 10]);
xlabel('SNR, E_s/N_0 (dB)'); ylabel('Bit error rate');
title('Uncoded BER: QPSK-GDSS (solid) vs QPSK-DSSS (dashed); markers = Monte Carlo');
legend(hG, arrayfun(@(n) sprintf('N = %d', n), Nlist, 'UniformOutput', false), 'Location', 'southwest');
save_fig(fh, cfg, 'fig09_ber_uncoded');
if cfg.saveFigs
    fid = fopen(fullfile(cfg.outdir, 'table_gdss_penalty.csv'), 'w');
    fprintf(fid, 'N,SNR_GDSS_1e-5_dB,SNR_DSSS_1e-5_dB,penalty_dB\n');
    fprintf(fid, '%d,%.2f,%.2f,%.2f\n', [Nlist; out.snrG; out.snrD; out.penalty]);
    fclose(fid);
end
fprintf('  Fig 9: GDSS penalty at BER 1e-5: %s dB (N = %s)\n', sprintf('%.2f ', out.penalty), sprintf('%d ', Nlist));
end
