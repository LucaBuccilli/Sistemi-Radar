"""Elaborazione a frame: Buffer -> finestra -> FFT -> spettro mediato.

È l'equivalente in Python della catena Simulink Buffer -> FFT -> |.|^2 -> media
(modello sim_fft_spectrum e blocco Spectrum Analyzer).
  - il picco cade nel bin k = f0 * N / Fs; la risoluzione è Fs / N
  - se f0 non cade su un bin, la finestra rettangolare "spande" energia (leakage), Hann molto meno
  - mediando M frame la varianza della stima di rumore cala di circa 1/M
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig, db

rng = np.random.default_rng(10)
Fs, N, M = 8000.0, 256, 64                     # campionamento, lunghezza frame, frame mediati
df = Fs / N
t = np.arange(N * M) / Fs
f_on = 1000.0                                  # esattamente sul bin 32
f_off = 1000.0 + df / 2                        # a metà tra due bin: caso peggiore di leakage
noise = 1e-4 * rng.standard_normal(t.size)       # rumore molto basso: rende visibile il leakage


def averaged_spectrum(x, win):
    frames = x.reshape(M, N)                   # Buffer: M frame da N campioni, senza sovrapposizione
    X = np.fft.rfft(frames * win, axis=1)
    return np.mean(np.abs(X) ** 2, axis=0) / np.sum(win**2)


f = np.fft.rfftfreq(N, 1 / Fs)
rect, hann = np.ones(N), np.hanning(N)
res = {}
for name, f0 in (("sul bin", f_on), ("tra due bin", f_off)):
    x = np.sin(2 * np.pi * f0 * t) + noise
    for wname, win in (("rettangolare", rect), ("Hann", hann)):
        S = averaged_spectrum(x, win)
        k = np.argmax(S)
        far = db(np.median(S[np.abs(f - f0) > 1000]) / S.max())      # livello lontano dal tono, relativo al picco
        res[(name, wname)] = S
        print(f"tono {name:11s} ({f0:7.1f} Hz), finestra {wname:12s}: picco al bin {k} ({f[k]:.1f} Hz), "
              f"livello a >1 kHz dal tono = {far:6.1f} dB")
print(f"Risoluzione in frequenza Fs/N = {df:.2f} Hz; bin atteso per 1000 Hz = {1000*N/Fs:.0f}")

# effetto della media sulla stima del rumore
x = noise
single = np.abs(np.fft.rfft(x[:N])) ** 2 / N
avg = averaged_spectrum(x, rect)
print(f"Rumore: deviazione standard relativa dello spettro, 1 frame = {np.std(single[5:-5])/np.mean(single[5:-5]):.2f}, "
      f"{M} frame = {np.std(avg[5:-5])/np.mean(avg[5:-5]):.3f} (atteso circa 1/sqrt(M) = {1/np.sqrt(M):.3f})")

fig, ax = plt.subplots(1, 2, figsize=(13, 4))
for (name, wname), S in res.items():
    if name == "tra due bin":
        ax[0].plot(f, db(S / S.max()), label=f"{wname}")
ax[0].set_xlim(0, 4000); ax[0].set_ylim(-110, 5); ax[0].grid(True); ax[0].legend()
ax[0].set_title(f"Tono a {f_off:.1f} Hz (tra due bin): leakage"); ax[0].set_xlabel("Hz"); ax[0].set_ylabel("dB")
ax[1].plot(f, db(single), alpha=0.6, label="1 frame")
ax[1].plot(f, db(avg), label=f"media di {M} frame")
ax[1].grid(True); ax[1].legend(); ax[1].set_xlabel("Hz"); ax[1].set_ylabel("dB")
ax[1].set_title("Solo rumore: la media riduce la varianza della stima")
savefig("d02_frame_fft_spectrum")
