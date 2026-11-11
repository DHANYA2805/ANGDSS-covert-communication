function out = fig01_dsss_theory(cfg)
%FIG01  Theoretical QPSK-DSSS BER for several spreading factors (paper Fig. 2).
snr_db = -50:0.25:20;
Nlist  = [1 32 256 1024];
fh = figure('Visible', cfg.visible, 'Position', [100 100 760 500]); hold on;
col = lines(numel(Nlist));
out.snr_at_1e5 = zeros(size(Nlist));
for i = 1:numel(Nlist)
    b = ber_dsss_theory(Nlist(i), 10.^(snr_db/10));
    semilogy(snr_db, max(b, 1e-300), 'LineWidth', 2, 'Color', col(i,:));
    out.snr_at_1e5(i) = required_snr(snr_db, b, 1e-5);
end
set(gca, 'YScale', 'log'); grid on; box on;
ylim([1e-7 1]); xlim([-50 20]);
xlabel('SNR, E_s/N_0 (dB)'); ylabel('Bit error rate');
title('Theoretical BER of QPSK-DSSS (Eq. 1)');
legend(arrayfun(@(n) sprintf('N = %d', n), Nlist, 'UniformOutput', false), 'Location', 'southwest');
save_fig(fh, cfg, 'fig01_dsss_theory_ber');
out.N = Nlist;
fprintf('  Fig 1: SNR for BER 1e-5: N=1 %.2f dB, N=256 %.2f dB (paper: 12.5 / -11.5 dB)\n', ...
    out.snr_at_1e5(1), out.snr_at_1e5(3));
end
