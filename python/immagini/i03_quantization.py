"""Quantizzazione: uniforme (regola dei 6 dB per bit) e ottima di Lloyd-Max.

Uniforme su [-A, A] con b bit: passo D = 2A/2^b, distorsione D^2/12 = sigma_x^2 / 2^(2b).
Lloyd-Max: soglie a metà tra i livelli, livelli = baricentri degli intervalli.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig, test_image

rng = np.random.default_rng(4)

# --- 1) quantizzatore uniforme
A = 1.0
x = rng.uniform(-A, A, 500_000)
print("bit | SNR misurato | teoria 6.02*b")
for b in range(2, 9):
    step = 2 * A / 2**b
    xq = step * (np.floor(x / step) + 0.5)                  # centro dell'intervallo
    snr = 10 * np.log10(np.var(x) / np.mean((x - xq) ** 2))
    print(f" {b}  |   {snr:6.2f} dB  |  {6.02*b:6.2f} dB")


# --- 2) Lloyd-Max per una gaussiana N(0,1)
def lloyd_max(samples, levels, iters=200):
    q = np.quantile(samples, (np.arange(levels) + 0.5) / levels)    # inizializzazione
    for _ in range(iters):
        theta = (q[:-1] + q[1:]) / 2                                # soglie a metà tra i livelli
        idx = np.searchsorted(theta, samples)
        q = np.array([samples[idx == k].mean() for k in range(levels)])   # baricentri
    theta = (q[:-1] + q[1:]) / 2
    dist = np.mean((samples - q[np.searchsorted(theta, samples)]) ** 2)
    return q, theta, dist


g = rng.standard_normal(1_000_000)
expected = {2: "+-0.7979", 3: "0, +-1.224", 4: "+-0.4528, +-1.5104"}
for L in (2, 3, 4):
    q, theta, dist = lloyd_max(g, L)
    print(f"Lloyd-Max {L} livelli: q = {np.round(q, 3)}, soglie = {np.round(theta, 3)}, "
          f"D = {dist:.4f}   (valori tabulati: {expected[L]})")

# confronto a 4 livelli: uniforme su [-3, 3] contro Lloyd-Max
step = 6 / 4
gu = np.clip(step * (np.floor(g / step) + 0.5), -2.25, 2.25)
print(f"Gaussiana, 4 livelli: D uniforme su [-3,3] = {np.mean((g-gu)**2):.4f}, D Lloyd-Max = {lloyd_max(g, 4)[2]:.4f}")

# --- 3) luminanza di un'immagine a 8, 4 e 2 bit
img = test_image(256)
fig, ax = plt.subplots(1, 3, figsize=(12, 4.2))
for a, b in zip(ax, (8, 4, 2)):
    L = 2**b
    a.imshow(np.floor(img * (L - 1) + 0.5) / (L - 1), cmap="gray", vmin=0, vmax=1)
    a.set_title(f"{b} bit ({L} livelli)" + (": falsi contorni" if b < 8 else "")); a.axis("off")
savefig("i03_quantization")
