function E = angdss_setup(cfg)
%ANGDSS_SETUP  Build the time-varying-noise test scenario used by Figs. 10-14.
%   E.g      : K x 1 noise multipliers (mean 1), slowly drifting log-normal
%   E.s2hat  : transmitter's noise estimates (median of L sensing samples)
%   E.sm2    : ANGDSS mask power per block (mean = GDSS power = 1)
%   E.sm2g   : ANGDSS mask power with perfect noise knowledge (genie)
%   E.tblG   : static-AWGN GDSS BER table for N (semi-analytical)
p   = cfg.env;
E.p = p;
E.g     = noise_environment(p.K, p.sigma_dB, p.a, p.clip_dB);
E.s2hat = estimate_noise_var(E.g, p.L, p.method, p.pImp, p.impGain_dB);
E.sm2   = angdss_mask_power(E.s2hat, 1, p.beta, p.peak);
E.sm2g  = angdss_mask_power(E.g, 1, p.beta, p.peak);
E.tblG.snr_db = -60:0.1:30;
E.tblG.ber    = ber_gdss_numeric(p.N, 10.^(E.tblG.snr_db / 10));
end
