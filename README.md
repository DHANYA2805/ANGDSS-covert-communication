# ANGDSS: Adaptive Noise-Matched Gaussian-Distributed Spread Spectrum

Covert (low-probability-of-detection) spread-spectrum communication in a drifting noise floor.

Gaussian-distributed spread spectrum (GDSS) masks every spread QPSK chip with the magnitude of a
Gaussian sample, so the transmitted samples look like thermal noise while the receiver still
despreads with the binary code alone. GDSS keeps the mask variance (the transmit power) fixed.
**ANGDSS** adds one control loop at the transmitter: before every block it senses the ambient
noise, estimates its power with a median estimator, and scales the mask variance in proportion,
under an average-power budget and an amplifier peak limit. The receiver is unchanged.

Base paper: I. Shakeel, J. Hilliard, W. Zhang, M. Rice, "Gaussian-Distributed Spread-Spectrum
for Covert Communications," *Sensors*, vol. 23, no. 8, 4081, 2023, doi:10.3390/s23084081.

## Key results (N = 256, noise floor drifting with 3 dB std, equal average power)

| Scheme | Average SNR for BER = 1e-5 |
|---|---|
| DSSS, fixed power | -6.95 dB |
| GDSS, fixed mask | -4.89 dB |
| **ANGDSS (median noise estimate)** | **-9.18 dB** |
| ANGDSS, exact noise knowledge | -9.33 dB |
| GDSS in constant noise (reference) | -9.34 dB |

- ANGDSS needs 4.3 dB less average power than fixed-mask GDSS for the same BER.
- At equal BER (1e-5), the warden can be 9.4 dB closer (in path loss) with ANGDSS than with GDSS
  for the same 95th-percentile detection probability of 0.2.
- The received amplitude distribution stays at the noise-only sampling floor (shape-based test).

## Repository layout

```
matlab/
  main_ANGDSS.m          settings + runs every figure
  functions/             transmitter, receiver, channel, BER theory, detectors,
                         noise sensing and mask-power controller
  figures/               one function per figure (fig01 ... fig14)
  results_expected/      figures, tables and summary.txt from the 'paper' setting
python/                  Python port of the functions, used for the extra analyses
                         (equal-BER detectability, causal controller,
                         feature-detector checks on ANGDSS)
  results/               figure produced by fig15.py
```

## How to run (MATLAB)

1. Open MATLAB in the `matlab/` folder.
2. In `main_ANGDSS.m` set `cfg.mode = 'quick'` (about 5-10 min) or `'paper'` (paper settings, 30-90 min).
3. Run. Figures, CSV tables, `results.mat` and `summary.txt` are written to `matlab/results/`.

Only base MATLAB (R2016b or newer) is needed; no toolboxes. The code also runs in GNU Octave 7+.
Random streams differ between MATLAB and Octave, so Monte Carlo values differ slightly from
`results_expected/summary.txt`; the semi-analytical curves are identical.

| Figure | Content |
|---|---|
| fig01-fig09 | Reproduction of the base paper: DSSS BER, constellations, Gaussianity and moments, 20th-moment detector, 4th-power detector, spectral coherence, number of convolutions, decision-variable PDF, uncoded BER |
| fig10 | ANGDSS loop: noise, estimate, mask power, chip SNR |
| fig11 | BER vs average SNR: DSSS, GDSS, ANGDSS (+ Monte Carlo) |
| fig12 | Energy-detector (radiometer) detectability |
| fig13 | Noise-matching (distribution-shape) score |
| fig14 | Gain vs drift depth; noise-estimator robustness |

## Python analyses

```
cd python
python3 verify.py    # rebuilds the drifting-noise scenario, checks the BER table, writes env.pkl
python3 extra.py     # equal-BER detectability and causal controller -> extra.json
python3 fig15.py     # equal-BER detectability + causal controller figure
python3 detchk.py    # 20th-moment test in a drifting floor
python3 cycchk.py    # cyclic-feature test on ANGDSS
```
Needs numpy, scipy, matplotlib.
