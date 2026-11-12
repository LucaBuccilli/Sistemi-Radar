"""DFT e convoluzione (circolare contro lineare, zero padding) e aliasing nel sottocampionamento.

1) Prodotto delle DFT = convoluzione CIRCOLARE. Con zero padding fino a M1+M2-1 si ottiene quella lineare.
2) Sottocampionare senza passa-basso produce aliasing; la decimazione (filtro + sottocampionamento) no.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage, signal
from utils import savefig, test_image

img = test_image(128)
h = np.ones((9, 9)) / 81.0                                  # filtro di media 9x9

# --- 1) convoluzione lineare, tramite DFT con e senza zero padding
lin = signal.convolve2d(img, h)                             # dimensione (128+9-1)^2
M = img.shape[0] + h.shape[0] - 1
dft_zp = np.real(np.fft.ifft2(np.fft.fft2(img, (M, M)) * np.fft.fft2(h, (M, M))))
circ = np.real(np.fft.ifft2(np.fft.fft2(img) * np.fft.fft2(h, img.shape)))
print(f"Dimensione uscita lineare: {lin.shape} = N1 + N2 - 1")
print(f"Errore massimo DFT con zero padding rispetto a conv2: {np.max(np.abs(lin - dft_zp)):.2e}")
print(f"Errore massimo DFT senza zero padding (circolare)   : "
      f"{np.max(np.abs(lin[:img.shape[0], :img.shape[1]] - circ)):.2f}  -> NON è la convoluzione lineare")

# --- 2) aliasing: zone plate (frequenza crescente dal centro)
n = 512
y, x = np.mgrid[-n // 2:n // 2, -n // 2:n // 2] / n
zone = 0.5 + 0.5 * np.cos(2 * np.pi * 120 * (x**2 + y**2))
Mfac = 4
down = zone[::Mfac, ::Mfac]                                 # solo sottocampionamento
dec = ndimage.gaussian_filter(zone, sigma=Mfac / 2)[::Mfac, ::Mfac]   # passa-basso, poi sottocampionamento

fig, ax = plt.subplots(1, 3, figsize=(12, 4.3))
for a, (im, t) in zip(ax, [(zone, "originale 512x512"),
                           (down, "sottocampionata x4: aliasing (falsi anelli)"),
                           (dec, "decimata x4 (passa-basso prima): nessun falso anello")]):
    a.imshow(im, cmap="gray", vmin=0, vmax=1); a.set_title(t, fontsize=9); a.axis("off")
savefig("i02_aliasing")

# --- 3) 1D: la DFT campiona la DTFT; un tono oltre Fs/2 si ripiega
Fs, f0, Nn = 8000.0, 5000.0, 800
xs = np.cos(2 * np.pi * f0 * np.arange(Nn) / Fs)
f = np.fft.rfftfreq(Nn, 1 / Fs)
print(f"Tono a {f0:.0f} Hz campionato a {Fs:.0f} Hz: picco dello spettro a "
      f"{f[np.argmax(np.abs(np.fft.rfft(xs)))]:.0f} Hz (atteso |f0 - Fs| = {abs(f0-Fs):.0f} Hz)")
