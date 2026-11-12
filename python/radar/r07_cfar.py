"""CA-CFAR (Cell Averaging CFAR) 2D sulla mappa range-Doppler.

Rivelatore quadratico: in assenza di bersaglio |z|^2 è esponenziale con media P_n (ignota e
variabile nella mappa). Il CA-CFAR stima P_n mediando N celle di addestramento attorno alla cella
sotto test (escludendo le celle di guardia) e usa la soglia
    T = alpha * P_n_stimata,    alpha = N * (Pfa^(-1/N) - 1)
che garantisce la Pfa desiderata qualunque sia il livello del disturbo.

1) Verifica della Pfa su rumore esponenziale i.i.d. con potenza variabile.
2) Catena completa su una mappa range-Doppler con clutter a Doppler zero e tre bersagli:
   soglia fissa, CA-CFAR 2D, CA-CFAR lungo la distanza e MTI + CA-CFAR.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
from utils import savefig, db

rng = np.random.default_rng(9)


def ca_cfar_2d(power, guard=(2, 2), train=(4, 8), pfa=1e-4):
    """CA-CFAR 2D. power: mappa |z|^2 (Doppler x distanza).
    guard/train: semiampiezza (Doppler, distanza) delle celle di guardia e di addestramento.
    Restituisce (rivelazioni, soglia, numero di celle di addestramento)."""
    gd, gr = guard
    td, tr = train
    kernel = np.ones((2 * (gd + td) + 1, 2 * (gr + tr) + 1))
    kernel[td:td + 2 * gd + 1, tr:tr + 2 * gr + 1] = 0           # toglie guardia + cella sotto test
    n_train = int(kernel.sum())
    alpha = n_train * (pfa ** (-1 / n_train) - 1)
    noise_est = convolve2d(power, kernel, mode="same", boundary="wrap") / n_train
    thr = alpha * noise_est
    return power > thr, thr, n_train


# --- 1) verifica della Pfa: rumore esponenziale con potenza che varia di 30 dB lungo la distanza
pfa = 1e-4
Nd, Nr = 256, 2048
level = 10 ** (np.linspace(0, 3, Nr) / 10)                         # da 0 a 30 dB
noise = rng.exponential(1.0, (Nd, Nr)) * level
det, _, n_train = ca_cfar_2d(noise, pfa=pfa)
alpha = n_train * (pfa ** (-1 / n_train) - 1)
fixed_thr = -np.log(pfa) * np.mean(level[: Nr // 10])             # soglia fissa tarata sulla zona a bassa potenza
print(f"Celle di addestramento N = {n_train}, alpha = {alpha:.2f} ({db(alpha):.1f} dB)")
print(f"Pfa progettata = {pfa:.1e}")
print(f"Pfa misurata CA-CFAR      = {det.mean():.2e}")
print(f"Pfa misurata soglia fissa = {(noise > fixed_thr).mean():.2e}  (il disturbo cresce: falsi allarmi)")
print(f"Perdita CFAR rispetto alla soglia ideale: {db(alpha / -np.log(pfa)):.2f} dB")


# --- 2) catena completa: compressione -> (MTI) -> FFT Doppler -> CA-CFAR
c, fc, B, tau, fs, PRT, Np = 3e8, 1.5e9, 3e6, 10e-6, 12e6, 100e-6, 32
lam = c / fc
t = np.arange(int(round(tau * fs))) / fs
chirp = np.exp(1j * np.pi * (B / tau) * (t - tau / 2) ** 2)
mf = chirp[::-1].conj()
ns = int(round(PRT * fs))
rng_axis = np.arange(ns) * c / (2 * fs)
vel_axis = np.fft.fftshift(np.fft.fftfreq(Np, PRT)) * lam / 2

data = (rng.standard_normal((Np, ns)) + 1j * rng.standard_normal((Np, ns))) / np.sqrt(2)
# clutter: echi fermi (Doppler 0) molto forti a breve distanza, decrescenti con la distanza
amp_cl = 30 * np.clip(1000 / np.maximum(rng_axis, 1), 0, 1) ** 1.5
cl_tx = amp_cl * (rng.standard_normal(ns) + 1j * rng.standard_normal(ns)) / np.sqrt(2)
data += np.convolve(cl_tx, chirp)[:ns][None, :]                    # uguale per tutti gli impulsi

targets = [(2500.0, 120.0, 0.15), (6000.0, 150.0, 0.5), (9000.0, -400.0, 0.1)]   # R, v, ampiezza
for R0, v, a in targets:
    d = int(round(2 * R0 / c * fs))
    fd = 2 * v / lam
    for n in range(Np):
        data[n, d:d + len(chirp)] += a * chirp * np.exp(1j * 2 * np.pi * fd * n * PRT)

rc = np.array([np.convolve(row, mf)[len(mf) - 1:] for row in data])   # compressione in distanza


def doppler_map(x):
    n = x.shape[0]
    return np.abs(np.fft.fftshift(np.fft.fft(x * np.hanning(n)[:, None], n=Np, axis=0), axes=0)) ** 2


P = doppler_map(rc)                                  # senza MTI
P_mti = doppler_map(rc[1:] - rc[:-1])                # canceller singolo prima della FFT

pfa2 = 1e-5
noise_floor = np.median(P[:, rng_axis > 12e3])
det_fixed = P > -np.log(pfa2) * noise_floor
det_2d, _, _ = ca_cfar_2d(P, guard=(2, 4), train=(4, 12), pfa=pfa2)        # finestra anche in Doppler
det_rng, _, _ = ca_cfar_2d(P, guard=(0, 6), train=(0, 16), pfa=pfa2)      # solo lungo la distanza
det_mti, thr, _ = ca_cfar_2d(P_mti, guard=(0, 6), train=(0, 16), pfa=pfa2)  # MTI + CFAR in distanza


def report(name, det):
    near = np.zeros_like(det)
    hits = []
    for R0, v, _ in targets:
        j = np.argmin(np.abs(rng_axis - R0)); i = np.argmin(np.abs(vel_axis - v))
        hits.append(bool(det[max(i - 1, 0):i + 2, max(j - 3, 0):j + 4].any()))
        near[max(i - 1, 0):i + 2, max(j - 10, 0):j + 11] = True
    ridge = np.abs(vel_axis)[:, None] < 40 * np.ones((1, ns))
    print(f"{name:22s}: bersagli {sum(hits)}/3 {hits}, celle sul clutter = {int((det & ridge).sum()):5d}, "
          f"altre false = {int((det & ~ridge & ~near).sum())}")


report("soglia fissa", det_fixed)
report("CA-CFAR 2D", det_2d)
report("CA-CFAR in distanza", det_rng)
report("MTI + CA-CFAR distanza", det_mti)

fig, ax = plt.subplots(2, 2, figsize=(13, 8), sharex=True, sharey=True)
ext = [rng_axis[0] / 1e3, rng_axis[-1] / 1e3, vel_axis[0], vel_axis[-1]]
ax[0, 0].imshow(db(P), aspect="auto", origin="lower", extent=ext, vmin=db(P).max() - 80)
ax[0, 0].set_title("Mappa range-Doppler [dB] (clutter a v = 0)")
ax[0, 1].imshow(det_fixed, aspect="auto", origin="lower", extent=ext, cmap="gray_r")
ax[0, 1].set_title(f"Soglia fissa: {det_fixed.sum()} celle")
ax[1, 0].imshow(det_rng, aspect="auto", origin="lower", extent=ext, cmap="gray_r")
ax[1, 0].set_title(f"CA-CFAR in distanza: {det_rng.sum()} celle")
ax[1, 1].imshow(det_mti, aspect="auto", origin="lower", extent=ext, cmap="gray_r")
ax[1, 1].set_title(f"MTI + CA-CFAR: {det_mti.sum()} celle")
for a in ax.ravel():
    for R0, v, _ in targets:
        a.plot(R0 / 1e3, v, "o", mfc="none", mec="r", ms=12)
for a in ax[1]:
    a.set_xlabel("distanza [km]")
for a in ax[:, 0]:
    a.set_ylabel("velocità [m/s]")
savefig("r07_cfar")
