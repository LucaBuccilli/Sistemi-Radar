"""Rivelazione a singolo impulso: soglia, Pfa, Pd (Monte Carlo) e formula di Albersheim.

H0: z = n            (inviluppo di Rayleigh)
H1: z = a + n        (inviluppo di Rice)
con n complesso gaussiano di potenza sigma_d^2 = 1 e SNR = |a|^2 / sigma_d^2.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig

rng = np.random.default_rng(0)


def noise(n):
    return (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)  # potenza 1


def albersheim(pfa, pd, n=1):
    """SNR richiesto [dB] per bersaglio non fluttuante, n impulsi integrati non coerentemente."""
    A = np.log(0.62 / pfa)
    Bq = np.log(pd / (1 - pd))
    return -5 * np.log10(n) + (6.2 + 4.54 / np.sqrt(n + 0.44)) * np.log10(A + 0.12 * A * Bq + 1.7 * Bq)


# 1) Verifica di Pfa = exp(-T^2 / sigma_d^2) con una Pfa "misurabile"
pfa = 1e-3
T = np.sqrt(-np.log(pfa))
pfa_mc = np.mean(np.abs(noise(2_000_000)) > T)
print(f"Pfa teorica = {pfa:.1e}, Pfa Monte Carlo = {pfa_mc:.2e}")

# 2) Curve Pd(SNR) per diverse Pfa
snr_db = np.arange(0, 18.5, 0.5)
plt.figure(figsize=(7, 4))
for pfa in (1e-2, 1e-4, 1e-6):
    T = np.sqrt(-np.log(pfa))
    pd = [np.mean(np.abs(np.sqrt(10 ** (s / 10)) + noise(200_000)) > T) for s in snr_db]
    plt.plot(snr_db, pd, label=f"Pfa = {pfa:.0e}")
    plt.plot(albersheim(pfa, 0.9), 0.9, "ko")
plt.xlabel("SNR [dB]"); plt.ylabel("Pd"); plt.grid(True); plt.legend()
plt.title("Pd in funzione dell'SNR (pallini: Albersheim per Pd = 0.9)")
savefig("r02_pd_vs_snr")

print(f"Albersheim: SNR per Pfa=1e-6, Pd=0.9, 1 impulso  = {albersheim(1e-6, 0.9):.1f} dB (atteso ~13 dB)")
print(f"Albersheim: SNR per Pfa=1e-6, Pd=0.9, 18 impulsi = {albersheim(1e-6, 0.9, 18):.1f} dB per impulso")
