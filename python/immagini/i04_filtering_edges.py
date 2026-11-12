"""Filtraggio spaziale 2D: smoothing, sharpening, gradiente (Sobel) e Laplaciano.

Verifica anche la risposta in frequenza della maschera a croce 3x3:
[0 1 0; 1 1 1; 0 1 0]  ->  H = 1 + 2 cos(w1) + 2 cos(w2).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from utils import savefig, test_image

rng = np.random.default_rng(5)
img = test_image(256)
noisy = img + 0.08 * rng.standard_normal(img.shape)

mean3 = np.ones((3, 3)) / 9
lap = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]], float)
sob_x = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], float)

smooth = ndimage.convolve(noisy, mean3, mode="nearest")           # passa-basso: toglie rumore
gauss = ndimage.gaussian_filter(noisy, 1.5)
laplacian = ndimage.convolve(img, lap, mode="nearest")            # derivata seconda: zero crossing sui bordi
sharp = img - 0.7 * laplacian                                     # originale + passa-alto
grad = np.hypot(ndimage.convolve(img, sob_x, mode="nearest"),
                ndimage.convolve(img, sob_x.T, mode="nearest"))   # modulo del gradiente: picco sui bordi

# Il passa-basso toglie rumore nelle zone lisce ma cancella anche i dettagli fini (le righe in basso):
# l'errore si misura quindi separatamente nelle due zone.
flat, tex = slice(0, 150), slice(170, 256)
for name, im in (("rumorosa", noisy), ("media 3x3", smooth), ("gaussiana sigma=1.5", gauss)):
    print(f"MSE {name:20s}: zona liscia = {np.mean((im[flat]-img[flat])**2):.5f}, "
          f"zona a righe = {np.mean((im[tex]-img[tex])**2):.5f}")

# risposta in frequenza della maschera a croce
cross = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], float)
w = np.linspace(-np.pi, np.pi, 65)
W1, W2 = np.meshgrid(w, w)
H = sum(cross[m + 1, n + 1] * np.exp(-1j * (W1 * n + W2 * m)) for m in (-1, 0, 1) for n in (-1, 0, 1))
print(f"Maschera a croce: errore rispetto a 1 + 2cos(w1) + 2cos(w2) = "
      f"{np.max(np.abs(H - (1 + 2*np.cos(W1) + 2*np.cos(W2)))):.1e}; H(0,0) = {H[32,32].real:.0f} (passa-basso)")

fig, ax = plt.subplots(2, 3, figsize=(12, 8))
items = [(noisy, "con rumore"), (smooth, "media 3x3 (smoothing)"), (sharp, "sharpening"),
         (grad, "gradiente (Sobel)"), (np.abs(laplacian), "|Laplaciano|"), (np.abs(H), "|H| maschera a croce")]
for a, (im, t) in zip(ax.ravel(), items):
    a.imshow(im, cmap="gray"); a.set_title(t); a.axis("off")
savefig("i04_filtering_edges")
