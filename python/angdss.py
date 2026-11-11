"""Python port of the ANGDSS MATLAB model (functions/ folder) for the extra
analyses requested in review: equal-BER detectability and a causal controller."""
import numpy as np
from scipy.special import erf, erfc, gammaincc, gammainc
from scipy.stats import gamma as gdist


def noise_environment(rng, K, sigma_dB, a=0.95, clip_dB=10):
    x = np.zeros(K)
    x[0] = sigma_dB * rng.standard_normal()
    e = rng.standard_normal(K)
    s = np.sqrt(1 - a * a) * sigma_dB
    for k in range(1, K):
        x[k] = a * x[k - 1] + s * e[k]
    x = np.clip(x, -clip_dB, clip_dB)
    g = 10 ** (x / 10)
    return g / g.mean()


def estimate_noise_var(rng, v, L, method='median', pImp=0.0, impGain_dB=20):
    K = v.size
    out = np.empty(K)
    B = max(1, 2 ** 21 // L)
    for k0 in range(0, K, B):
        idx = slice(k0, min(K, k0 + B))
        n = v[idx].size
        wI = rng.standard_normal((L, n)) / np.sqrt(2)
        wQ = rng.standard_normal((L, n)) / np.sqrt(2)
        if pImp > 0:
            hit = rng.random((L, n)) < pImp
            a = np.sqrt(10 ** (impGain_dB / 10) / 2)
            wI = wI + hit * a * rng.standard_normal((L, n))
            wQ = wQ + hit * a * rng.standard_normal((L, n))
        e = (wI ** 2 + wQ ** 2) * v[idx][None, :]
        out[idx] = e.mean(0) if method == 'mean' else np.median(e, 0) / np.log(2)
    return out


def mask_power(s2hat, Pbar=1.0, beta=1.0, peak=4.0):
    rho = Pbar / s2hat.mean()
    sm2 = beta * rho * s2hat
    for _ in range(20):
        sm2 = np.minimum(sm2, peak * Pbar)
        sm2 = sm2 * (beta * Pbar / sm2.mean())
        if sm2.max() <= peak * Pbar * (1 + 1e-9):
            break
    return sm2


def mask_power_causal(s2hat, Pbar=1.0, peak=4.0, tau=200, warm=None):
    """Causal controller: reference noise level is an exponentially weighted
    average of PAST estimates only; no renormalisation over future blocks."""
    lam = 1.0 / tau
    m = s2hat[0] if warm is None else warm
    sm2 = np.empty_like(s2hat)
    for k in range(s2hat.size):
        sm2[k] = min(Pbar * s2hat[k] / m, peak * Pbar)
        m = (1 - lam) * m + lam * s2hat[k]
    return sm2


def ber_dsss(N, snr):
    return 0.5 * erfc(np.sqrt(N * snr / 2))


def pdf_chip(x, sn, ss=1.0):
    v = ss ** 2 + sn ** 2
    return np.exp(-x ** 2 / (2 * v)) * (1 + erf(ss * x / (sn * np.sqrt(2 * v)))) / np.sqrt(2 * np.pi * v)


def ber_gdss_numeric(N, snr_lin):
    snr_lin = np.atleast_1d(snr_lin)
    ber = np.zeros(snr_lin.size)
    bits = [int(b) for b in bin(N)[2:][::-1]]
    for k, s in enumerate(snr_lin):
        sn = np.sqrt(1 / s)
        mu1 = np.sqrt(2 / np.pi)
        v1 = (1 - 2 / np.pi) + sn ** 2
        muZ, sdZ = N * mu1, np.sqrt(N * v1)
        sr = np.sqrt(1 + sn ** 2)
        dx0 = min(sr, sn * sr) / 16
        W = 28 * sdZ
        Lg = int(2 ** np.ceil(np.log2(np.ceil(W / dx0))))
        Lg = min(max(Lg, 2048), 2 ** 21)
        dx = W / Lg
        x1 = mu1 + dx * (np.arange(Lg) - Lg / 2)
        y1 = pdf_chip(x1, sn) * dx
        y1 /= y1.sum()
        P = np.fft.fft(y1)
        acc = None
        for b in bits:
            if b:
                acc = P.copy() if acc is None else acc * P
            P = P * P
        yN = np.real(np.fft.ifft(acc))
        yN[yN < 0] = 0
        yN /= yN.sum()
        z0 = N * x1[0] + dx * np.arange(Lg)
        zz = z0 - W * np.floor((z0 - (muZ - W / 2)) / W)
        ber[k] = yN[zz < 0].sum()
    return ber


class Table:
    def __init__(self, N, lo=-60, hi=30, step=0.1):
        self.snr_db = np.arange(lo, hi + 1e-9, step)
        self.lb = np.log10(np.maximum(ber_gdss_numeric(N, 10 ** (self.snr_db / 10)), 1e-300))

    def __call__(self, snr_lin):
        s = 10 * np.log10(np.maximum(snr_lin, 1e-30))
        s = np.clip(s, self.snr_db[0], self.snr_db[-1])
        return 10 ** np.interp(s, self.snr_db, self.lb)


def _hist(x):
    q = np.round(x * 100).astype(int)
    u, inv = np.unique(q, return_inverse=True)
    w = np.bincount(inv) / x.size
    return u / 100, w


def ber_curves(sa_db, g, sm2, N, tbl):
    cF, wF = _hist(-10 * np.log10(g))
    cA, wA = _hist(10 * np.log10(sm2 / g))
    out = {k: np.zeros(sa_db.size) for k in ('dsss', 'gdss', 'angdss', 'adsss')}
    for i, s in enumerate(sa_db):
        gF = 10 ** ((s + cF) / 10)
        gA = 10 ** ((s + cA) / 10)
        out['dsss'][i] = wF @ ber_dsss(N, gF)
        out['gdss'][i] = wF @ tbl(gF)
        out['angdss'][i] = wA @ tbl(gA)
        out['adsss'][i] = wA @ ber_dsss(N, gA)
    return out


def required_snr(sa, ber, t=1e-5):
    lb = np.log10(np.maximum(ber, 1e-300))
    lt = np.log10(t)
    i = np.where((lb[:-1] >= lt) & (lb[1:] < lt))[0]
    if i.size == 0:
        return np.nan
    i = i[0]
    return sa[i] + (lt - lb[i]) * (sa[i + 1] - sa[i]) / (lb[i + 1] - lb[i])


def pd_radiometer(gam, nw, pfa=0.05):
    tau = gdist.ppf(1 - pfa, nw)
    return gdist.sf(tau / (1 + gam), nw)
