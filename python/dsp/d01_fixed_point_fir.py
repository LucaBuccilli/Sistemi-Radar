"""Filtro FIR in virgola fissa: effetto della lunghezza di parola su uscita e risposta in frequenza.

Formato Q con segno: w bit totali, f bit frazionari -> passo 2^-f, range [-2^(w-f-1), 2^(w-f-1)).
  ingresso e uscita : f = w - 2   (range +-2: somma di due toni di ampiezza 1)
  coefficienti      : f = w - 1   (|h| < 1)
  accumulatore      : precisione piena (come nei DSP, che hanno accumulatori più larghi)
Stesso filtro del modello Simulink sim_fir_fixed_point: passa-basso a 1 kHz, Fs = 8 kHz, 65 coefficienti.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz, lfilter
from utils import savefig, db


def quantize(x, w, f):
    """Arrotonda al passo 2^-f e satura nel range del formato con w bit (con segno)."""
    q = 2.0 ** -f
    lo, hi = -(2.0 ** (w - f - 1)), 2.0 ** (w - f - 1) - q
    return np.clip(np.round(x / q) * q, lo, hi)


Fs, fc, Nt = 8000.0, 1000.0, 64
n = np.arange(Nt + 1) - Nt / 2
h = 2 * fc / Fs * np.sinc(2 * fc / Fs * n) * np.hamming(Nt + 1)   # sinc finestrata

t = np.arange(8000) / Fs
x = np.sin(2 * np.pi * 300 * t) + np.sin(2 * np.pi * 2500 * t)
y_ref = lfilter(h, 1, x)                                            # riferimento in doppia precisione

print("bit | SNR uscita [dB] | attenuazione in banda oscura [dB]")
results = {}
for w in (8, 10, 12, 14, 16, 20):
    hq = quantize(h, w, w - 1)
    yq = quantize(lfilter(hq, 1, quantize(x, w, w - 2)), w, w - 2)
    snr = db(np.mean(y_ref**2) / np.mean((y_ref - yq) ** 2))
    f_ax, H = freqz(hq, 1, 4096, fs=Fs)
    stop = -20 * np.log10(np.max(np.abs(H[f_ax > 1500])))
    results[w] = (snr, stop, f_ax, H)
    print(f" {w:2d} |     {snr:6.1f}      |     {stop:6.1f}")
s8, s16 = results[8][0], results[16][0]
print(f"Pendenza misurata: {(s16 - s8) / 8:.2f} dB per bit (teoria circa 6.02 dB/bit)")
_, H_ref = freqz(h, 1, 4096, fs=Fs)
print(f"Attenuazione in banda oscura in doppia precisione: {-20*np.log10(np.max(np.abs(H_ref[f_ax > 1500]))):.1f} dB")

fig, ax = plt.subplots(1, 2, figsize=(13, 4))
for w in (8, 12, 16):
    ax[0].plot(results[w][2], 20 * np.log10(np.abs(results[w][3]) + 1e-12), label=f"{w} bit")
ax[0].plot(f_ax, 20 * np.log10(np.abs(H_ref) + 1e-12), "k--", label="double")
ax[0].set_ylim(-110, 5); ax[0].set_xlabel("Hz"); ax[0].set_ylabel("|H| [dB]"); ax[0].grid(True); ax[0].legend()
ax[0].set_title("Coefficienti quantizzati: la banda oscura risale")
ws = sorted(results)
ax[1].plot(ws, [results[w][0] for w in ws], "o-", label="misurato")
ax[1].plot(ws, [s8 + 6.02 * (w - 8) for w in ws], "k--", label="6.02 dB/bit")
ax[1].set_xlabel("lunghezza di parola [bit]"); ax[1].set_ylabel("SNR uscita [dB]"); ax[1].grid(True); ax[1].legend()
ax[1].set_title("SNR del FIR in virgola fissa")
savefig("d01_fixed_point_fir")
