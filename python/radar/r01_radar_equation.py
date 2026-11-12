"""Equazione radar: SNR in funzione della distanza, portata massima, orizzonte radar.

Caso di riferimento: radar di sorveglianza in banda S (lambda = 0.2 m) e dimensionamento della potenza.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig, db

k, T0, c = 1.380649e-23, 290.0, 3e8

# --- Caso di riferimento: Pt=32 kW, G=35 dB, lambda=0.2 m, RCS=10 m^2, L=10 dB, F=5 dB, B=1 MHz
Pt, G, lam, sigma = 32e3, 10 ** 3.5, 0.2, 10.0
L, F, B = 10.0, 10 ** 0.5, 1e6


def snr(R):
    """SNR = Pt G^2 lambda^2 sigma / ((4 pi)^3 R^4 L k T0 B F)."""
    return Pt * G**2 * lam**2 * sigma / ((4 * np.pi) ** 3 * R**4 * L * k * T0 * B * F)


for R in (100, 1e3, 10e3, 100e3, 400e3):
    print(f"R = {R/1e3:7.1f} km   SNR = {db(snr(R)):6.1f} dB")

snr_min = 10 ** (13.3 / 10)  # Pd = 0.9, Pfa = 1e-6, bersaglio non fluttuante
Rmax = (Pt * G**2 * lam**2 * sigma / ((4 * np.pi) ** 3 * L * k * T0 * B * F * snr_min)) ** 0.25
print(f"Portata massima con SNR_min = 13.3 dB: {Rmax/1e3:.1f} km")

R = np.logspace(2, np.log10(400e3), 400)
plt.figure(figsize=(7, 4))
plt.semilogx(R / 1e3, db(snr(R)), label="SNR singolo impulso")
plt.axhline(13.3, color="r", ls="--", label="SNR richiesto (13.3 dB)")
plt.axvline(Rmax / 1e3, color="gray", ls=":")
plt.xlabel("Distanza [km]"); plt.ylabel("SNR [dB]"); plt.grid(True, which="both"); plt.legend()
plt.title("Equazione radar: l'SNR cala di 40 dB per decade (1/R^4)")
savefig("r01_snr_vs_range")

# --- Dimensionamento: potenza di picco per Rmax = 150 km
f, Bx, eta, Ag, Fx, sig, Rx, snr_req = 2e9, 1e6, 0.6, 12.0, 10 ** 0.4, 5.0, 150e3, 10 ** 1.2
lam_x = c / f
Gx = 4 * np.pi * eta * Ag / lam_x**2
Pt_req = snr_req * (4 * np.pi) ** 3 * Rx**4 * k * T0 * Fx * Bx / (Gx**2 * lam_x**2 * sig)
print(f"Dimensionamento: G = {db(Gx):.2f} dB, Pt richiesta = {Pt_req/1e3:.1f} kW (calcolo a mano in dB: 87.7 kW)")

# --- Orizzonte radar con raggio terrestre equivalente 4/3
Re = 4 / 3 * 6371e3
for hA in (100.0, 4.0):
    row = [np.sqrt(2 * Re * hA) + np.sqrt(2 * Re * hT) for hT in (10e3, 900.0, 100.0, 0.0)]
    print(f"h_antenna = {hA:5.0f} m -> orizzonte [km] per h_T = 10 km, 900 m, 100 m, 0: "
          + ", ".join(f"{r/1e3:.0f}" for r in row))
