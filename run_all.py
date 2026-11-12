"""Lancia in sequenza tutte le simulazioni Python e salva le figure in python/figures/.

Uso (dalla cartella della repository):
    python run_all.py            tutte le simulazioni
    python run_all.py r04 i05    solo quelle il cui nome inizia con r04 o i05
"""
import glob
import os
import subprocess
import sys
import time

root = os.path.dirname(os.path.abspath(__file__))
scripts = sorted(glob.glob(os.path.join(root, "python", "radar", "r*.py"))) + \
          sorted(glob.glob(os.path.join(root, "python", "immagini", "i*.py")))
filters = sys.argv[1:]
if filters:
    scripts = [s for s in scripts if any(os.path.basename(s).startswith(f) for f in filters)]

failed = []
t0 = time.time()
for s in scripts:
    print("\n" + "=" * 70 + f"\n{os.path.basename(s)}\n" + "=" * 70, flush=True)
    if subprocess.run([sys.executable, s], cwd=root).returncode != 0:
        failed.append(os.path.basename(s))

print(f"\nCompletate {len(scripts) - len(failed)}/{len(scripts)} simulazioni in {time.time() - t0:.0f} s.")
print(f"Figure in: {os.path.join(root, 'python', 'figures')}")
if failed:
    print("Errori in: " + ", ".join(failed))
    sys.exit(1)
