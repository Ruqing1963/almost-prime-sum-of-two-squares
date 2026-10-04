import numpy as np
_p = __file__.replace("richert_fine.py", "richert.py")
_src = open(_p, encoding="utf-8-sig").read().split("deltas = np.arange(0.08")[0]
_src = _src.replace('__file__.replace("richert.py", "net_yield.py")', repr(_p.replace("richert.py", "net_yield.py")))
_src = _src.replace('__file__.replace("richert.py", "sieve_funcs.py")', repr(_p.replace("richert.py", "sieve_funcs.py")))
exec(_src)

for K in (4.0, 5.0):
    b = scan(K, True, np.arange(0.04, 0.16, 0.01), np.arange(0.002, 0.03, 0.001))
    for k in sorted(b)[:3]:
        d, a, lam, b1, net = b[k]
        print(f"K={K}: Richert k={k}: delta={d:.2f}, alpha={a:.3f}, lam={lam:.3f}, beta1={b1:.2f}, net kept={net:+.4f}")
