"""Trasformata di Fourier 2D: spettro di un'immagine, passa-basso e passa-alto ideali, rotazione.

Proprietà verificata: ruotare l'immagine ruota lo spettro dello stesso angolo.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from utils import savefig, test_image

img = test_image(256)
N = img.shape[0]


def spectrum(a):
    """Modulo dello spettro in scala logaritmica, con la continua al centro."""
    return np.log1p(np.abs(np.fft.fftshift(np.fft.fft2(a))))


F = np.fft.fftshift(np.fft.fft2(img))
u = np.arange(-N // 2, N // 2)
U, V = np.meshgrid(u, u)
radius = np.hypot(U, V)
lp = np.real(np.fft.ifft2(np.fft.ifftshift(F * (radius < 20))))    # passa-basso ideale
hp = np.real(np.fft.ifft2(np.fft.ifftshift(F * (radius >= 20))))   # passa-alto ideale
rot = ndimage.rotate(img, 30, reshape=False, mode="nearest")

# Le righe oblique dell'immagine danno due picchi in (+-24, +-24) cicli per immagine
Fa = np.abs(F).copy(); Fa[N // 2, N // 2] = 0
Fa[radius < 5] = 0
iy, ix = np.unravel_index(np.argmax(Fa), Fa.shape)
print(f"Picco dello spettro (esclusa la continua): ({u[ix]}, {u[iy]}) cicli/immagine; atteso (+-24, +-24)")
print(f"Energia nelle basse frequenze (raggio < 20): {np.sum(np.abs(F[radius<20])**2)/np.sum(np.abs(F)**2)*100:.1f} %")

fig, ax = plt.subplots(2, 4, figsize=(14, 7))
for a, (im, title) in zip(ax[0], [(img, "originale"), (lp, "passa-basso (sfocata, con ringing)"),
                                  (hp, "passa-alto (solo bordi)"), (rot, "ruotata di 30 gradi")]):
    a.imshow(im, cmap="gray"); a.set_title(title); a.axis("off")
for a, im in zip(ax[1], [img, lp, hp, rot]):
    a.imshow(spectrum(im), cmap="gray"); a.set_title("spettro (log)"); a.axis("off")
savefig("i01_fourier2d")
