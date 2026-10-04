"""
Figures for the repository (illustrative, floating point; not part of the certified computation).
  figures/semilinear_sieve_functions.png : F_{1/2}, f_{1/2} (kappa = 1/2, beta = 1)
  figures/smallest_k_by_scheme.png       : smallest k found by the floating-point scans for each scheme
Run from the repository root:  python code/make_figures.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXP = os.path.join(HERE, "exploration")

src = open(os.path.join(EXP, "richert.py"), encoding="utf-8-sig").read().split("deltas = np.arange(0.08")[0]
src = src.replace('__file__.replace("richert.py", "net_yield.py")', repr(os.path.join(EXP, "net_yield.py")))
src = src.replace('__file__.replace("richert.py", "sieve_funcs.py")', repr(os.path.join(EXP, "sieve_funcs.py")))
ns = {"__file__": os.path.join(EXP, "richert.py")}
exec(src, ns)

os.makedirs(os.path.join(ROOT, "figures"), exist_ok=True)

s = np.linspace(0.3, 5, 600)
fig, ax = plt.subplots(figsize=(6, 3.6))
ax.plot(s, ns["Fs"](s), label=r"$F_{1/2}(s)$")
ax.plot(s, ns["fs"](s), label=r"$f_{1/2}(s)$")
ax.axhline(1, color="0.6", lw=0.8)
ax.axvline(1, color="0.6", lw=0.8, ls=":")
ax.set_xlabel("s"); ax.set_ylim(0, 2.8); ax.legend(frameon=False)
ax.set_title(r"Semi-linear sieve functions ($\kappa=1/2$, $\beta=1$): $f_{1/2}(s)=0$ for $s\leq 1$", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "semilinear_sieve_functions.png"), dpi=200)

labels, values = [], []
for K in (4.0, 4.5, 5.0):
    best = ns["scan"](K, False, np.arange(0.08, 0.24, 0.01), np.arange(0.010, 0.12, 0.0005))
    labels.append(f"no weights\nK={K:g}"); values.append(min(best))
for K in (4.0, 5.0):
    best = ns["scan"](K, True, np.arange(0.04, 0.16, 0.01), np.arange(0.002, 0.03, 0.001))
    labels.append(f"Richert\nK={K:g}"); values.append(min(best))
fig, ax = plt.subplots(figsize=(6, 3.6))
bars = ax.bar(labels, values, color=["0.65"] * 3 + ["C0"] * 2)
ax.bar_label(bars)
ax.set_ylim(0, max(values) + 5)
ax.set_ylabel(r"smallest $k$ with positive net yield")
ax.set_title("Floating-point parameter scans (illustrative; only K=4, Richert, k=9 is certified)", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "figures", "smallest_k_by_scheme.png"), dpi=200)
print("smallest k:", dict(zip([l.replace(chr(10), ' ') for l in labels], values)))
print("figures written to", os.path.join(ROOT, "figures"))
