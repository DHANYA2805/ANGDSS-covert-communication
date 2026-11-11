function out = fig12_angdss_covertness(cfg, E)
%FIG12  Warden detectability (energy detector + KL divergence) over time.
%   Worst case for the transmitter: the warden sees the same channel SNR as
%   the receiver and knows the current noise power exactly (ideal
%   radiometer). Observation window n_w chips, false-alarm rate P_FA.
%   T = sum|y|^2 / sigma_n^2 ~ Gamma(n_w,1) under H0 and (1+g)Gamma(n_w,1)
%   under H1, so P_D = 1 - gammainc(tau/(1+g), n_w).
%   KL divergence per window D = n_w (g - ln(1+g)); any detector obeys
%   P_FA + P_MD >= 1 - sqrt(D/2)  (Bash et al., covert-communication bound).
nw  = cfg.env.nw;  PFA = cfg.env.PFA;
tau = gammaincinv(1 - PFA, nw);
pd  = @(gam) 1 - gammainc(tau ./ (1 + gam), nw);
kl  = @(gam) nw * (gam - log(1 + gam));
sa0 = cfg.env.demoSNR;  nv0 = 10^(-sa0/10);
kk  = 1:min(400, numel(E.g));
gG  = 1 ./ (nv0 * E.g);  gA = E.sm2 ./ (nv0 * E.g);
pdG = pd(gG);  pdA = pd(gA);
% Monte Carlo check of the radiometer at the demo SNR
nRep = cfg.env.edReps; hitG = 0; hitA = 0; nTot = 0;
for k = kk
    wG = cgauss(nw * nRep, E.g(k) * nv0);
    yG = gdss_mask(dsss_spread(qpsk_mod(nw*nRep/cfg.env.N), 2*randi([0 1], cfg.env.N, 1)-1), 1) + wG;
    yA = sqrt(E.sm2(k)) * gdss_mask(dsss_spread(qpsk_mod(nw*nRep/cfg.env.N), 2*randi([0 1], cfg.env.N, 1)-1), 1) ...
         + cgauss(nw * nRep, E.g(k) * nv0);
    TG = sum(reshape(abs(yG).^2, nw, nRep), 1) / (E.g(k) * nv0);
    TA = sum(reshape(abs(yA).^2, nw, nRep), 1) / (E.g(k) * nv0);
    hitG = hitG + sum(TG > tau);  hitA = hitA + sum(TA > tau);  nTot = nTot + nRep;
end
out.mcPdGDSS = hitG / nTot;  out.mcPdANGDSS = hitA / nTot;
out.thPdGDSS = mean(pdG(kk)); out.thPdANGDSS = mean(pdA(kk));
% versus average SNR
sa = -30:0.5:0;
mG = zeros(size(sa)); mA = mG; qG = mG; qA = mG; dG = mG; dA = mG;
for i = 1:numel(sa)
    nv = 10^(-sa(i)/10);
    a = pd(1 ./ (nv * E.g));  b = pd(E.sm2 ./ (nv * E.g));
    mG(i) = mean(a);  mA(i) = mean(b);
    sG = sort(a); sA = sort(b); j = ceil(0.95 * numel(sG));
    qG(i) = sG(j);    qA(i) = sA(j);
    dG(i) = max(kl(1 ./ (nv * E.g)));  dA(i) = max(kl(E.sm2 ./ (nv * E.g)));
end
fh = figure('Visible', cfg.visible, 'Position', [60 60 1250 470]);
subplot(1,2,1); hold on;
plot(kk, pdG(kk), '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.4);
plot(kk, pdA(kk), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
plot(kk, PFA * ones(size(kk)), 'k:');
grid on; box on; xlim([kk(1) kk(end)]); ylim([0 1]);
xlabel('Block index k'); ylabel('Detection probability P_D');
title(sprintf('(a) Ideal radiometer, n_w = %d, P_{FA} = %.2f, avg SNR = %d dB', nw, PFA, sa0));
legend('GDSS (fixed mask)', 'ANGDSS', 'P_{FA}', 'Location', 'northwest');
subplot(1,2,2); hold on;
plot(sa, mG, '-',  'Color', [0.00 0.45 0.74], 'LineWidth', 2);
plot(sa, qG, '--', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.6);
plot(sa, mA, '-',  'Color', [0.85 0.33 0.10], 'LineWidth', 2);
plot(sa, qA, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
grid on; box on; ylim([0 1]); xlim([sa(1) sa(end)]);
xlabel('Average SNR (dB)'); ylabel('Detection probability P_D');
title('(b) Mean and 95th-percentile P_D over all blocks');
legend('GDSS mean', 'GDSS 95th pct', 'ANGDSS mean', 'ANGDSS 95th pct', 'Location', 'northwest');
save_fig(fh, cfg, 'fig12_angdss_covertness');
% covertness budget: highest average SNR with 95th-percentile P_D <= 0.2
lim = 0.2;
out.saCovertGDSS   = interp1(qG + 1e-12*(1:numel(qG)), sa, lim);
out.saCovertANGDSS = interp1(qA + 1e-12*(1:numel(qA)), sa, lim);
out.sa = sa; out.meanG = mG; out.meanA = mA; out.q95G = qG; out.q95A = qA;
out.maxKLG = dG; out.maxKLA = dA;
out.allPdGDSS = mean(pdG); out.allPdANGDSS = mean(pdA);
fprintf('  Fig 12: avg SNR %d dB, all blocks: mean P_D GDSS %.3f, ANGDSS %.3f\n', sa0, out.allPdGDSS, out.allPdANGDSS);
fprintf('          first %d blocks: theory GDSS %.3f / MC %.3f, ANGDSS %.3f / MC %.3f\n', ...
    numel(kk), out.thPdGDSS, out.mcPdGDSS, out.thPdANGDSS, out.mcPdANGDSS);
fprintf('          highest avg SNR with 95%% of blocks at P_D <= %.1f: GDSS %.2f dB, ANGDSS %.2f dB\n', ...
    lim, out.saCovertGDSS, out.saCovertANGDSS);
end
