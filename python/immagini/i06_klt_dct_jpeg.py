"""KLT, DCT e schema JPEG (DCT a blocchi 8x8 + quantizzazione con la matrice di luminanza).

1) KLT: autovettori della covarianza, coefficienti scorrelati (esempio 2x2 risolvibile a mano).
2) DCT: base ortonormale, approssima la KLT dei processi di Markov molto correlati.
3) JPEG: blocchi 8x8 -> DCT -> round(X / Q) -> la maggior parte dei coefficienti va a zero.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig, test_image, psnr

rng = np.random.default_rng(7)

# --- 1) KLT, esempio 2x2: Kx = [2 1; 1 2]
Kx = np.array([[2.0, 1.0], [1.0, 2.0]])
lam, U = np.linalg.eigh(Kx)
lam, U = lam[::-1], U[:, ::-1]                    # ordine decrescente
xs = rng.multivariate_normal([0, 0], Kx, 200_000)
Xk = xs @ U                                       # coefficienti KLT
print(f"KLT: autovalori = {lam} (attesi: 3 e 1), primo autovettore = {np.round(np.abs(U[:,0]), 4)}")
print(f"     covarianza dei coefficienti =\n{np.round(np.cov(Xk.T), 3)}  (diagonale: scorrelati)")
print(f"     energia nel primo coefficiente: {lam[0]/lam.sum()*100:.0f} %")


# --- 2) matrice DCT NxN: righe = funzioni di base
def dct_matrix(n=8):
    k, m = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    Cm = np.sqrt(2 / n) * np.cos(np.pi * (2 * m + 1) * k / (2 * n))
    Cm[0, :] = 1 / np.sqrt(n)
    return Cm


Cm = dct_matrix(8)
print(f"DCT: errore di ortonormalità = {np.max(np.abs(Cm @ Cm.T - np.eye(8))):.1e}")

# DCT come KLT di un processo di Markov con rho = 0.95
rho = 0.95
Km = rho ** np.abs(np.subtract.outer(np.arange(8), np.arange(8)))
_, Uk = np.linalg.eigh(Km)
Uk = Uk[:, ::-1]
print("Somiglianza |<base DCT k, autovettore KLT k>| per k = 0..7 (Markov, rho = 0.95):",
      np.round(np.abs(np.sum(Cm * Uk.T, axis=1)), 4))

# --- 3) JPEG semplificato
Q = np.array([[16, 11, 10, 16, 24, 40, 51, 61], [12, 12, 14, 19, 26, 58, 60, 55],
              [14, 13, 16, 24, 40, 57, 69, 56], [14, 17, 22, 29, 51, 87, 80, 62],
              [18, 22, 37, 56, 68, 109, 103, 77], [24, 35, 55, 64, 81, 104, 113, 92],
              [49, 64, 78, 87, 103, 121, 120, 101], [72, 92, 95, 98, 112, 100, 103, 99]], float)
img = test_image(256) * 255


def jpeg_like(image, scale):
    out = np.zeros_like(image)
    nonzero = 0
    for r in range(0, image.shape[0], 8):
        for c in range(0, image.shape[1], 8):
            X = Cm @ (image[r:r + 8, c:c + 8] - 128) @ Cm.T      # DCT 2D separabile
            Xq = np.round(X / (Q * scale))                       # quantizzazione (qui si perde informazione)
            nonzero += np.count_nonzero(Xq)
            out[r:r + 8, c:c + 8] = Cm.T @ (Xq * Q * scale) @ Cm + 128
    return out, nonzero / image.size


fig, ax = plt.subplots(1, 4, figsize=(15, 4))
ax[0].imshow(img, cmap="gray", vmin=0, vmax=255); ax[0].set_title("originale"); ax[0].axis("off")
for a, s in zip(ax[1:], (1, 4, 12)):
    rec, frac = jpeg_like(img, s)
    p = psnr(img / 255, rec / 255)
    print(f"JPEG-like, Q x {s:2d}: coefficienti non nulli = {frac*100:5.1f} %, PSNR = {p:.1f} dB")
    a.imshow(rec, cmap="gray", vmin=0, vmax=255)
    a.set_title(f"Q x {s}: {frac*100:.0f}% coeff., {p:.0f} dB"); a.axis("off")
savefig("i06_jpeg")

# le 64 funzioni di base della DCT 8x8
fig, ax = plt.subplots(8, 8, figsize=(6, 6))
for k1 in range(8):
    for k2 in range(8):
        ax[k1, k2].imshow(np.outer(Cm[k1], Cm[k2]), cmap="gray"); ax[k1, k2].axis("off")
savefig("i06_dct_basis")
