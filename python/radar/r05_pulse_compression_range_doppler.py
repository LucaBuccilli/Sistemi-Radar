"""Compressione d'impulso (chirp + filtro adattato) e mappa range-Doppler di un treno coerente.

Schema: filtro adattato in tempo veloce (distanza), poi FFT in tempo lento (Doppler).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig, db

c = 3e8
fc, B, tau, fs = 1.5e9, 3e6, 10e-6, 12e6          # prodotto tempo-banda B*tau = 30
PRT, N = 100e-6, 32
lam = c / fc
rng = np.random.default_rng(3)

t = np.arange(int(round(tau * fs))) / fs
chirp = np.exp(1j * np.pi * (B / tau) * (t - tau / 2) ** 2)     # LFM in banda base
mf = chirp[::-1].conj()                                         # filtro adattato h(t) = s*(-t)

# 1) Compressione del singolo impulso
comp = np.abs(np.convolve(chirp, mf))
comp /= comp.max()
width = np.sum(comp > 1 / np.sqrt(2)) / fs
print(f"Durata impulso: {tau*1e6:.1f} us -> dopo la compressione (-3 dB): {width*1e6:.2f} us "
      f"(teoria 0.886/B = {0.886e6/B:.2f} us, misura limitata dal passo 1/fs = {1e6/fs:.2f} us)")
print(f"Risoluzione in distanza: {c*tau/2:.0f} m senza compressione, {c/(2*B):.0f} m con compressione")

# 2) Treno di N impulsi: due bersagli + rumore
ns = int(round(PRT * fs))
targets = [(6000.0, 150.0, 1.0), (9000.0, -400.0, 0.5)]        # distanza [m], velocità [m/s], ampiezza
data = (rng.standard_normal((N, ns)) + 1j * rng.standard_normal((N, ns))) / np.sqrt(2)
for R0, v, a in targets:
    d = int(round(2 * R0 / c * fs))
    fd = 2 * v / lam
    for n in range(N):
        data[n, d:d + len(chirp)] += a * chirp * np.exp(1j * 2 * np.pi * fd * n * PRT)

rc = np.array([np.convolve(row, mf)[len(mf) - 1:] for row in data])   # compressione in distanza
rd = np.fft.fftshift(np.fft.fft(rc * np.hanning(N)[:, None], axis=0), axes=0)  # FFT Doppler

rng_axis = np.arange(ns) * c / (2 * fs)
vel_axis = np.fft.fftshift(np.fft.fftfreq(N, PRT)) * lam / 2
power = db(np.abs(rd) ** 2)
i, j = np.unravel_index(np.argmax(power), power.shape)
print(f"Picco della mappa: R = {rng_axis[j]:.0f} m, v = {vel_axis[i]:.0f} m/s (bersaglio vero: 6000 m, 150 m/s; "
      f"passo dei filtri Doppler {lam/(2*N*PRT):.0f} m/s)")
print(f"Velocità non ambigua: +/- {lam/(4*PRT):.0f} m/s, risoluzione in velocità: {lam/(2*N*PRT):.0f} m/s")

plt.figure(figsize=(8, 4))
plt.pcolormesh(rng_axis / 1e3, vel_axis, power, shading="auto", vmin=power.max() - 50)
plt.colorbar(label="dB"); plt.xlabel("Distanza [km]"); plt.ylabel("Velocità radiale [m/s]")
plt.title("Mappa range-Doppler (filtro adattato + FFT in tempo lento)")
savefig("r05_range_doppler")

plt.figure(figsize=(7, 3.5))
lag = (np.arange(len(comp)) - (len(chirp) - 1)) / fs * 1e6
plt.plot(lag, 20 * np.log10(np.maximum(comp, 1e-4)))
plt.xlabel("Ritardo [us]"); plt.ylabel("dB"); plt.grid(True); plt.ylim(-50, 1)
plt.title("Uscita del filtro adattato al chirp (lobi laterali a circa -13 dB)")
savefig("r05_pulse_compression")
