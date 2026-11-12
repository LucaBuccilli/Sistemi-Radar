"""Integrazione di N impulsi: coerente (FFT) contro non coerente (somma dei moduli quadri).

I rivelatori derivano dal test GLRT (rapporto di verosimiglianza generalizzato):
  coerente:      | sum_n r_n exp(-j 2 pi fd n T) |   -> guadagno N
  non coerente:  sum_n |r_n|^2                      -> guadagno minore di N
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig

rng = np.random.default_rng(1)
N, trials, pfa = 16, 40_000, 1e-3
kd = 3                                           # bin Doppler del bersaglio: fd = kd / (N T)
steer = np.exp(1j * 2 * np.pi * kd * np.arange(N) / N)


def noise():
    return (rng.standard_normal((trials, N)) + 1j * rng.standard_normal((trials, N))) / np.sqrt(2)


stats = {
    "singolo impulso": lambda r: np.abs(r[:, 0]),
    "non coerente (N=16)": lambda r: np.sum(np.abs(r) ** 2, axis=1),
    "coerente (N=16)": lambda r: np.abs(r @ steer.conj()),
}
# soglie fissate sul solo rumore, in modo da avere la stessa Pfa per tutti
thr = {name: np.quantile(f(noise()), 1 - pfa) for name, f in stats.items()}

snr_db = np.arange(-12, 16, 1.0)
pd = {name: [] for name in stats}
for s in snr_db:
    r = np.sqrt(10 ** (s / 10)) * steer + noise()
    for name, f in stats.items():
        pd[name].append(np.mean(f(r) > thr[name]))

plt.figure(figsize=(7, 4))
for name in stats:
    plt.plot(snr_db, pd[name], label=name)
plt.axhline(0.9, color="gray", ls=":")
plt.xlabel("SNR per impulso [dB]"); plt.ylabel("Pd"); plt.grid(True); plt.legend()
plt.title(f"Integrazione di impulsi (Pfa = {pfa:.0e})")
savefig("r03_pulse_integration")

need = {name: float(np.interp(0.9, pd[name], snr_db)) for name in stats}
for name, v in need.items():
    print(f"SNR per impulso per Pd=0.9, {name:22s}: {v:5.1f} dB")
print(f"Guadagno coerente     : {need['singolo impulso'] - need['coerente (N=16)']:.1f} dB "
      f"(teorico 10log10(N) = {10*np.log10(N):.1f} dB)")
print(f"Guadagno non coerente : {need['singolo impulso'] - need['non coerente (N=16)']:.1f} dB")
