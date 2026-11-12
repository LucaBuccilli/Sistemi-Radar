"""Image recovery: y = h * x + n. Filtro inverso, pseudo-inverso, minimi quadrati vincolati e Wiener.

Tutti i filtri sono applicati in frequenza (la sfocatura è una convoluzione circolare):
  inverso          F = 1 / H
  pseudo-inverso   F = 1 / H dove |H| > soglia, 0 altrove
  CLS              F = H* / (|H|^2 + lambda |C|^2),  C = Laplaciano (passa-alto)
  Wiener           F = H* / (|H|^2 + NSR),           NSR = Pww / Pxx
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from utils import savefig, test_image, psnr

rng = np.random.default_rng(6)
x = test_image(256)
N = x.shape[0]

# PSF: sfocatura da movimento lungo la diagonale (15 pixel), centrata nell'origine
psf = np.zeros((N, N))
L = 15
for i in range(L):
    psf[(i - L // 2) % N, (i - L // 2) % N] = 1 / L
H = np.fft.fft2(psf)

sigma_n = 0.01
y = np.real(np.fft.ifft2(np.fft.fft2(x) * H)) + sigma_n * rng.standard_normal(x.shape)
Y = np.fft.fft2(y)

lap = np.zeros((N, N))
lap[0, 0], lap[0, 1], lap[1, 0], lap[0, -1], lap[-1, 0] = 4, -1, -1, -1, -1
C = np.fft.fft2(lap)

Pxx = np.abs(np.fft.fft2(x)) ** 2 / N**2          # spettro di potenza del segnale (qui noto: caso ideale)
Pww = sigma_n**2 * np.ones((N, N))

Hs = np.where(np.abs(H) < 1e-12, 1e-12, H)        # evita la divisione per zero nel filtro inverso
filters = {
    "inverso": 1 / Hs,
    "pseudo-inverso": np.where(np.abs(H) > 0.1, 1 / Hs, 0),
    "CLS (lambda = 0.01)": H.conj() / (np.abs(H) ** 2 + 0.01 * np.abs(C) ** 2),
    "Wiener (NSR costante)": H.conj() / (np.abs(H) ** 2 + sigma_n**2 / np.var(x)),
    "Wiener (NSR vero)": H.conj() / (np.abs(H) ** 2 + Pww / np.maximum(Pxx, 1e-12)),
}
results = {name: np.real(np.fft.ifft2(Y * F)) for name, F in filters.items()}

print(f"PSNR immagine osservata      : {psnr(x, y):6.2f} dB")
for name, r in results.items():
    print(f"PSNR {name:24s}: {psnr(x, r):6.2f} dB")

fig, ax = plt.subplots(2, 3, figsize=(12, 8))
items = [(x, "originale"), (y, "sfocata + rumore")] + [(results[k], k) for k in
         ("inverso", "pseudo-inverso", "CLS (lambda = 0.01)", "Wiener (NSR vero)")]
for a, (im, t) in zip(ax.ravel(), items):
    a.imshow(np.clip(im, 0, 1), cmap="gray", vmin=0, vmax=1); a.set_title(t); a.axis("off")
savefig("i05_restoration")
