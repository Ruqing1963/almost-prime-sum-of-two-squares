"""
Round 7: robustness of the switching constant K, and Richert logarithmic weights.

Weighted sum over survivors n (class-3 primes < y = N^(1/2-delta) and class-1 primes < z0 = N^alpha sifted):
   w_n = 1 - lam * sum_{q | n, q = 1 mod 4, z0 <= q < y1} (1 - log q / log y1),   y1 = N^beta1
For a good survivor (no class-3 factor) with w_n > 0:  Omega(n) < 1/lam + 1/beta1  (and always < 1/alpha).
Lower bound (normalised by Sing(N) N / log^2 N):
   e^-gamma (L - lam R) / (2 sqrt(alpha (1/2-delta)))  -  K c(delta, alpha)
   L = best vector-sieve lower function at level 1/2
   R = (1/2) int_alpha^beta1 (1 - u/beta1) U(u) du/u,   U(u) = min over level split of F3 F1 at level 1/2 - u
       (upper vector sieve a b <= a+ b+ applied to A_q, q = N^u; class-1 primes have density 1/2)
The bad survivors (two class-3 factors >= y) have w_n <= 1, so subtracting K c (their switching bound) is valid.
"""
import numpy as np
_src = open(__file__.replace("richert.py", "net_yield.py"), encoding="utf-8-sig").read().split("print(\"sanity")[0]
_src = _src.replace('__file__.replace("net_yield.py", "sieve_funcs.py")', repr(__file__.replace("richert.py", "sieve_funcs.py")))
exec(_src)

A12 = 2 * np.sqrt(EG / np.pi)


def Fs(s):
    s = np.asarray(s, dtype=float)
    return np.where(s <= 2, A12 / np.sqrt(np.maximum(s, 1e-12)), at(Fh, np.minimum(s, S_MAX)))


def fs(s):
    return at(fh, np.minimum(np.asarray(s, dtype=float), S_MAX))


def L_best(delta, alpha):
    t3 = np.linspace(0.001, 0.499, 499)
    s3, s1 = t3 / (0.5 - delta), (0.5 - t3) / alpha
    return float(np.max(fs(s3) * Fs(s1) + Fs(s3) * fs(s1) - Fs(s3) * Fs(s1)))


def U_curve(delta, alpha, u):
    out = np.empty_like(u)
    for i, ui in enumerate(u):
        t3 = np.linspace(0.0005, 0.5 - ui - 0.0005, 300)
        out[i] = np.min(Fs(t3 / (0.5 - delta)) * Fs((0.5 - ui - t3) / alpha))
    return out


def scan(K, richert, deltas, alphas, betas=np.arange(0.10, 0.50, 0.01)):
    best = {}
    for delta in deltas:
        for alpha in alphas:
            L = L_best(delta, alpha)
            norm = 2 * np.sqrt(alpha * (0.5 - delta)) * np.exp(np.euler_gamma)
            slack = L - K * c_switch(delta, alpha) * norm      # = L - lam R must stay above this
            if slack <= 0:
                continue
            cands = [(int(np.ceil(1 / alpha - 1e-9)) - 1, 0.0, None, slack / norm)]  # lam = 0
            if richert:
                u = np.linspace(alpha, 0.495, 160)
                U = U_curve(delta, alpha, u)
                for b1 in betas:
                    if b1 <= alpha:
                        continue
                    m = u <= b1
                    R = 0.5 * np.trapz((1 - u[m] / b1) * U[m] / u[m], u[m])
                    lam = 0.95 * slack / R                       # keep 5% of the slack as margin
                    bound = min(1 / alpha, 1 / lam + 1 / b1)
                    cands.append((int(np.ceil(bound - 1e-9)) - 1, lam, b1, 0.05 * slack / norm))
            for k, lam, b1, net in cands:
                if k not in best or net > best[k][-1]:
                    best[k] = (delta, alpha, lam, b1, net)
    return best


deltas = np.arange(0.08, 0.24, 0.01)
for K in (4.0, 4.5, 5.0):
    b = scan(K, False, deltas, np.arange(0.010, 0.12, 0.0005))
    k0 = min(b)
    d, a, _, _, net = b[k0]
    print(f"K={K}: no weights   -> smallest k = {k0}  (delta={d:.2f}, alpha={a:.4f}, net={net:+.4f})"
          f";  net at k=28: {b.get(28, (0,0,0,0,float('nan')))[-1]:+.4f}, k=30: {b.get(30, (0,0,0,0,float('nan')))[-1]:+.4f}")

for K in (4.0, 5.0):
    b = scan(K, True, np.arange(0.08, 0.24, 0.02), np.arange(0.01, 0.12, 0.0025))
    k0 = min(b)
    d, a, lam, b1, net = b[k0]
    print(f"K={K}: Richert      -> smallest k = {k0}  (delta={d:.2f}, alpha={a:.4f}, lam={lam:.3f}, beta1={b1}, "
          f"net kept={net:+.4f})")
