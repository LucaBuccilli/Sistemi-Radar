"""Canceller MTI binomiali: risposta in frequenza, attenuazione del clutter e Improvement Factor.

Clutter con spettro gaussiano: rho(kT) = exp(-2 pi^2 sigma_f^2 k^2 T^2).
Canceller di ordine M: H(f) = (1 - exp(-j 2 pi f T))^M, |H|^2 = (4 sin^2(pi f T))^M.
Scenario: radar di ricerca a corto raggio in banda S, sigma_v = 0.2 m/s, lambda = 0.1 m, PRF = 1250 Hz.
"""
import os, sys
from math import comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter
from utils import savefig, db

PRF, lam, sigma_v = 1250.0, 0.1, 0.2
T = 1 / PRF
sigma_f = 2 * sigma_v / lam                       # 4 Hz


def rho(k, sf=sigma_f):
    return np.exp(-2 * np.pi**2 * sf**2 * (k * T) ** 2)


def clutter_out(M, sf=sigma_f):
    """Potenza di clutter in uscita (ingresso unitario): sum_m (-1)^m C(2M, M+m) rho(mT)."""
    return sum((-1) ** m * comb(2 * M, M + m) * rho(m, sf) for m in range(-M, M + 1))


def improvement_factor(fd, M, sf=sigma_f):
    return (4 * np.sin(np.pi * fd * T) ** 2) ** M / clutter_out(M, sf)


print(f"sigma_f = {sigma_f:.1f} Hz, rho(T) = {rho(1):.4f}, rho(2T) = {rho(2):.4f}, rho(3T) = {rho(3):.4f}")
print("Canceller | CA [dB] | IF a PRF/2 | IF a PRF/20 | IF medio")
for M in (1, 2, 3):
    ca = db(1 / clutter_out(M))
    if_avg = db(comb(2 * M, M) / clutter_out(M))   # guadagno medio sul bersaglio = C(2M, M)
    print(f"  M = {M}   | {ca:6.1f}  | {db(improvement_factor(PRF/2, M)):8.1f}   |"
          f" {db(improvement_factor(PRF/20, M)):8.1f}    | {if_avg:6.1f}")

# Effetto della rotazione dell'antenna: lo spettro del clutter si allarga.
# Nota: arrotondando rho(T) a 0.999 (come spesso si fa nei conti a mano) si otterrebbero
# circa 87.5 dB; con il valore esatto (0.9996) il risultato è più alto.
ToT = 66.7e-3
sigma_rot = np.hypot(sigma_f, 0.265 / ToT)
print(f"Con rotazione dell'antenna: sigma_f' = {sigma_rot:.2f} Hz -> IF triplo canceller a PRF/2 = "
      f"{db(improvement_factor(PRF/2, 3, sigma_rot)):.1f} dB")

fd = np.linspace(0, PRF, 1000)
plt.figure(figsize=(7, 4))
for M in (1, 2, 3):
    plt.plot(fd, db(improvement_factor(fd, M)), label=f"ordine {M}")
plt.ylim(-20, 120); plt.xlabel("Doppler del bersaglio [Hz]"); plt.ylabel("Improvement Factor [dB]")
plt.title("IF dei canceller binomiali: velocità cieche a 0 e PRF"); plt.grid(True); plt.legend()
savefig("r04_mti_improvement_factor")

# Verifica per simulazione: clutter gaussiano correlato + canceller singolo
rng = np.random.default_rng(2)
n = 2**18
f = np.fft.fftfreq(n, T)
shape = np.exp(-f**2 / (4 * sigma_f**2))          # radice della PSD gaussiana
cl = np.fft.ifft(shape * (rng.standard_normal(n) + 1j * rng.standard_normal(n)))
cl /= np.sqrt(np.mean(np.abs(cl) ** 2))
out = lfilter([1, -1], 1, cl)[1:]                 # y[n] = x[n] - x[n-1]
print(f"Simulazione canceller singolo: CA = {db(1/np.mean(np.abs(out)**2)):.1f} dB "
      f"(teoria 1/(2(1-rho)) = {db(1/clutter_out(1)):.1f} dB)")
