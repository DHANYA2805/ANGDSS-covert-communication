%% =====================================================================
%  ANGDSS: Adaptive Noise-Matched Gaussian-Distributed Spread Spectrum
%  ---------------------------------------------------------------------
%  Part A  - Reproduction of the base paper
%            I. Shakeel, J. Hilliard, W. Zhang, M. Rice, "Gaussian-Distributed
%            Spread-Spectrum for Covert Communications", Sensors 23, 4081, 2023
%            Fig. 1  DSSS theoretical BER                  (paper Fig. 2)
%            Fig. 2  DSSS / GDSS constellations            (paper Figs. 3, 5)
%            Fig. 3  GDSS distribution + moments table     (paper Fig. 7, Table 1)
%            Fig. 4  20th-order-moment detector            (paper Fig. 8)
%            Fig. 5  Modulation-stripping detector         (paper Fig. 9)
%            Fig. 6  Spectral-coherence detector           (paper Fig. 10)
%            Fig. 7  Number of numerical convolutions      (paper Fig. 12)
%            Fig. 8  Decision-variable PDF                 (paper Fig. 13)
%            Fig. 9  Uncoded BER, GDSS vs DSSS             (paper Fig. 14)
%  Part B  - Proposed ANGDSS (our novelty) in a drifting noise floor
%            Fig. 10 Adaptive loop: sensing, estimate, mask power, chip SNR
%            Fig. 11 BER vs average SNR: DSSS / GDSS / ANGDSS (+ Monte Carlo)
%            Fig. 12 Warden detectability (radiometer P_D, KL divergence)
%            Fig. 13 Noise-matching stealth score (KL to ambient-noise PDF)
%            Fig. 14 Gain vs fluctuation depth, estimator robustness
%
%  Needs only base MATLAB (R2016b or newer). No Communications or
%  Statistics Toolbox is required. Also runs in GNU Octave >= 7.
%  Run time (typical laptop): 'quick' ~5-10 min, 'paper' ~30-90 min
%  (the 20th-moment test with 1e6 bits x 5 trials dominates).
%  Outputs: PNG (and .fig) figures, CSV tables, results.mat and
%  summary.txt in the folder ./results
%% =====================================================================
clc; clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(fullfile(here, 'functions'), fullfile(here, 'figures'));

%% ----------------------- user settings -------------------------------
cfg.mode     = 'quick';          % 'quick' (fast check) or 'paper' (full settings)
cfg.seed     = 2026;             % random seed -> reproducible results
cfg.outdir   = fullfile(here, 'results');
cfg.saveFigs = true;
cfg.visible  = 'on';             % 'off' to run without showing figures
runPartA     = true;
runPartB     = true;

% ANGDSS scenario (Part B)
cfg.env.N          = 256;        % spreading factor
cfg.env.Ns         = 32;         % QPSK symbols per block (one mask update per block)
cfg.env.K          = 100000;     % blocks in the noise trace
cfg.env.sigma_dB   = 3;          % std of the noise-floor drift (dB)
cfg.env.a          = 0.95;       % block-to-block correlation of the drift
cfg.env.clip_dB    = 10;         % drift limited to +/- 10 dB
cfg.env.L          = 256;        % noise samples sensed per block
cfg.env.method     = 'median';   % 'median' (robust), 'mean' or 'genie'
cfg.env.pImp       = 0.01;       % impulsive interference in sensing slot
cfg.env.impGain_dB = 20;         % impulse power above the noise floor
cfg.env.beta       = 1;          % safety knob of sigma_m^2 = beta*rho*sigma_n^2
cfg.env.peak       = 4;          % HPA limit: sigma_m^2 <= 4 x average (+6 dB)
cfg.env.demoSNR    = -15;        % average SNR for the time-trace figures
cfg.env.nw         = 1024;       % warden observation window (chips)
cfg.env.PFA        = 0.05;       % warden false-alarm probability
cfg.env.sigmaSweep = 0:1:6;      % Fig. 14a
cfg.env.Lsweep     = [8 16 32 64 128 256 512 1024];   % Fig. 14b

switch cfg.mode
    case 'paper'
        cfg.tableBits    = 1e6;  cfg.momentBits = 1e6;  cfg.momentTrials = 5;
        cfg.momentSNR    = -20:2:0;
        cfg.modstripSym  = 10000; cfg.scfSym = 4096;    cfg.pdfTrials = 1e5;
        cfg.mcN = [16 64 256];   cfg.mcMinBER = 1e-5;  cfg.mcMaxBits = 1e7;  cfg.mcMinErr = 100;
        cfg.env.mcSNR = -20:2:-10; cfg.env.mcMaxBits = 1e7; cfg.env.edReps = 200;
        cfg.env.Kw = 4000; cfg.env.nbW = 512;
    otherwise   % quick
        cfg.tableBits    = 2e5;  cfg.momentBits = 1e6;  cfg.momentTrials = 2;
        cfg.momentSNR    = -20:4:0;
        cfg.modstripSym  = 10000; cfg.scfSym = 2048;    cfg.pdfTrials = 5e4;
        cfg.mcN = [16 64 256];   cfg.mcMinBER = 1e-4;  cfg.mcMaxBits = 2e6;  cfg.mcMinErr = 100;
        cfg.env.mcSNR = -20:3:-11; cfg.env.mcMaxBits = 4e6; cfg.env.edReps = 50;
        cfg.env.Kw = 2000; cfg.env.nbW = 512;
end

%% ----------------------- run -----------------------------------------
rng(cfg.seed);
set(0, 'DefaultAxesFontSize', 11, 'DefaultLineLineWidth', 1.2);
if ~exist(cfg.outdir, 'dir'), mkdir(cfg.outdir); end
diary(fullfile(cfg.outdir, 'summary.txt')); diary on;
fprintf('ANGDSS simulation, mode = %s, seed = %d, %s\n', cfg.mode, cfg.seed, datestr(now));
R = struct();
t0 = tic;
if runPartA
    fprintf('\n===== Part A: base paper (GDSS vs DSSS) =====\n');
    rng(cfg.seed + 1); R.f01 = fig01_dsss_theory(cfg);
    rng(cfg.seed + 2); R.f02 = fig02_constellations(cfg);
    rng(cfg.seed + 3); R.f03 = fig03_gaussianity(cfg);
    rng(cfg.seed + 4); R.f04 = fig04_moment_detector(cfg);
    rng(cfg.seed + 5); R.f05 = fig05_mod_stripping(cfg);
    rng(cfg.seed + 6); R.f06 = fig06_spectral_correlation(cfg);
    rng(cfg.seed + 7); R.f07 = fig07_convolutions(cfg);
    rng(cfg.seed + 8); R.f08 = fig08_decision_pdf(cfg);
    rng(cfg.seed + 9); R.f09 = fig09_ber_uncoded(cfg);
    fprintf('  Part A time: %.1f s\n', toc(t0));
end
if runPartB
    fprintf('\n===== Part B: proposed ANGDSS =====\n');
    rng(cfg.seed + 100); E = angdss_setup(cfg);     % same scenario in every mode
    rng(cfg.seed + 10); R.f10 = fig10_angdss_tracking(cfg, E);
    rng(cfg.seed + 11); R.f11 = fig11_angdss_ber(cfg, E);
    rng(cfg.seed + 12); R.f12 = fig12_angdss_covertness(cfg, E);
    rng(cfg.seed + 13); R.f13 = fig13_angdss_noise_match(cfg, E);
    rng(cfg.seed + 14); R.f14 = fig14_angdss_sensitivity(cfg, E);
end
fprintf('\nTotal run time: %.1f s\n', toc(t0));
save(fullfile(cfg.outdir, 'results.mat'), 'R', 'cfg', '-v7');
diary off;
