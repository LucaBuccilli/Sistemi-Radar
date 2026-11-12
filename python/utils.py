"""Funzioni comuni a tutti gli script Python della repository."""
import os

import matplotlib
matplotlib.use("Agg")  # nessuna finestra: le figure vengono salvate su file
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")


def savefig(name):
    """Salva la figura corrente in python/figures/<name>.png."""
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, name + ".png")
    plt.tight_layout()
    plt.savefig(path, dpi=120)
    plt.close()
    print(f"[figura salvata] {os.path.relpath(path)}")


def db(x):
    """Potenza lineare -> dB."""
    return 10 * np.log10(np.maximum(x, 1e-300))


def test_image(n=256):
    """Immagine di prova sintetica in [0,1]: sfumatura, rettangolo, disco e righe oblique.

    Si usa un'immagine sintetica per non dipendere da file esterni.
    """
    y, x = np.mgrid[0:n, 0:n] / n
    img = 0.3 + 0.3 * x                                   # sfumatura (basse frequenze)
    img[(x > 0.15) & (x < 0.45) & (y > 0.2) & (y < 0.5)] = 0.9   # rettangolo (bordi netti)
    img[(x - 0.7) ** 2 + (y - 0.35) ** 2 < 0.02] = 0.1          # disco
    stripes = 0.5 + 0.4 * np.sin(2 * np.pi * 24 * (x + y))      # tessitura a righe oblique
    mask = y > 0.65
    img[mask] = stripes[mask]
    return img


def psnr(ref, x):
    """PSNR in dB per immagini in [0,1]."""
    return 10 * np.log10(1.0 / np.mean((ref - x) ** 2))
