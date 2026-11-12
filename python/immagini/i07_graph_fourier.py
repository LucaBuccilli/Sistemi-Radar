"""Segnali su grafo: Laplaciano L = D - A e trasformata di Fourier su grafo (GFT).

Base = autovettori di L; gli autovalori fanno la parte delle frequenze.
Un segnale "liscio" sul grafo ha i coefficienti concentrati sugli autovalori piccoli.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig

rng = np.random.default_rng(8)

# grafo con 3 gruppi (cluster) di 10 nodi: molti archi dentro i gruppi, pochi tra gruppi
centers = np.array([[0, 0], [3, 0.5], [1.5, 2.5]])
pos = np.vstack([c + 0.5 * rng.standard_normal((10, 2)) for c in centers])
label = np.repeat(np.arange(3), 10)
n = len(pos)
A = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        p = 0.7 if label[i] == label[j] else 0.02
        A[i, j] = A[j, i] = rng.random() < p
for a, b in ((0, 10), (10, 20), (20, 0)):          # garantisce che il grafo sia connesso
    A[a, b] = A[b, a] = 1

D = np.diag(A.sum(axis=1))
Lap = D - A
lam, U = np.linalg.eigh(Lap)                       # autovalori crescenti = "frequenze" crescenti
print(f"Autovalori più piccoli: {np.round(np.abs(lam[:4]), 3)} (il primo è 0: autovettore costante)")

smooth = np.array([1.0, 3.0, -2.0])[label] + 0.1 * rng.standard_normal(n)   # quasi costante nei gruppi
rough = rng.standard_normal(n)                                            # varia a caso tra nodi vicini
for name, s in (("liscio", smooth), ("irregolare", rough)):
    S = U.T @ s                                    # GFT: S_k = <s, u_k>
    e3 = np.sum(S[:3] ** 2) / np.sum(S**2) * 100
    print(f"Segnale {name:10s}: s'Ls = {s @ Lap @ s:7.1f}, energia nei primi 3 coefficienti = {e3:5.1f} %")
print(f"Ricostruzione s = U S: errore = {np.max(np.abs(U @ (U.T @ smooth) - smooth)):.1e}")

fig, ax = plt.subplots(1, 3, figsize=(14, 4))
for i in range(n):
    for j in range(i + 1, n):
        if A[i, j]:
            ax[0].plot(pos[[i, j], 0], pos[[i, j], 1], color="lightgray", lw=0.6, zorder=1)
sc = ax[0].scatter(pos[:, 0], pos[:, 1], c=smooth, cmap="jet", zorder=2)
plt.colorbar(sc, ax=ax[0]); ax[0].set_title("Segnale liscio sul grafo"); ax[0].axis("off")
ax[1].stem(lam, np.abs(U.T @ smooth)); ax[1].set_title("GFT del segnale liscio"); ax[1].set_xlabel("autovalore")
ax[2].stem(lam, np.abs(U.T @ rough)); ax[2].set_title("GFT del segnale irregolare"); ax[2].set_xlabel("autovalore")
savefig("i07_graph_fourier")
