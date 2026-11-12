"""SAR: chirp in azimut e focalizzazione. La risoluzione lungo traccia vale circa d_a / 2.

Scenario avionico in banda X: lambda = 3.1 cm, antenna d_a = 1.8 m, R = 20 km.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig

lam, da, Ry, v = 0.031, 1.8, 20e3, 100.0
Ds = lam * Ry / da                      # impronta in azimut (apertura sintetica)
Tobs = Ds / v                           # tempo di osservazione di un punto
Bd = 2 * v / da                         # banda Doppler
PRF = 1.5 * Bd                          # sopra il minimo 2v/d_a
print(f"Impronta reale: {Ds:.0f} m, Tobs = {Tobs:.2f} s, banda Doppler = {Bd:.1f} Hz, PRF = {PRF:.1f} Hz")
print(f"Pendenza del chirp: {2*v**2/(lam*Ry):.2f} Hz/s, prodotto tempo-banda = {Bd*Tobs:.0f}")

ta = np.arange(-2 * Tobs, 2 * Tobs, 1 / PRF)     # tempo lento
x = v * ta                                       # posizione della piattaforma


def echo(x0):
    """Eco di un punto in x0: fase -4 pi R(t)/lambda, visibile solo dentro il fascio."""
    R = np.sqrt(Ry**2 + (x - x0) ** 2)
    return np.exp(-1j * 4 * np.pi / lam * (R - Ry)) * (np.abs(x - x0) <= Ds / 2)


ref = echo(0.0)
sig = echo(0.0) + echo(2.0)                      # due punti distanti 2 m lungo traccia


def focus(s, x0):
    """Filtro adattato al punto in x0: correlazione con l'eco atteso (compressione in azimut)."""
    return np.abs(np.vdot(echo(x0), s))


# Il passo tra impulsi è v/PRF = 0.6 m: per misurare la risoluzione si valuta
# l'uscita su una griglia più fitta di posizioni x0.
xf = np.arange(-6, 8, 0.02)
foc = np.array([focus(sig, x0) for x0 in xf])
foc /= foc.max()
single = np.array([focus(ref, x0) for x0 in xf])
single /= single.max()
res = np.ptp(xf[single > 1 / np.sqrt(2)])
print(f"Risoluzione reale (senza SAR): {Ds:.0f} m")
print(f"Risoluzione focalizzata misurata (-3 dB): {res:.2f} m; teoria 0.886 * d_a/2 = {0.886*da/2:.2f} m "
      f"(ordine di grandezza d_a/2 = {da/2:.2f} m)")

fig, ax = plt.subplots(1, 3, figsize=(12, 3.5))
ax[0].plot(ta, np.real(ref)); ax[0].set_title("Eco di un punto: chirp in tempo lento")
ax[0].set_xlabel("tempo lento [s]"); ax[0].set_xlim(-Tobs / 2, Tobs / 2)
fi = np.gradient(np.unwrap(np.angle(ref))) * PRF / (2 * np.pi)
m = np.abs(ta) < Tobs / 2
ax[1].plot(ta[m], fi[m]); ax[1].set_title("Frequenza Doppler istantanea (lineare)")
ax[1].set_xlabel("tempo lento [s]"); ax[1].set_ylabel("Hz"); ax[1].grid(True)
ax[2].plot(xf, foc); ax[2].set_xlim(-6, 8); ax[2].grid(True)
ax[2].set_title("Dopo la focalizzazione: due punti a 2 m"); ax[2].set_xlabel("lungo traccia [m]")
savefig("r06_sar_azimuth")
