"""
Net yield for  N = p + n,  n = 1 mod 4 a sum of two squares, all prime factors >= z0 = N^alpha  (=> Omega(n) < 1/alpha).

Sieve on A = {N - p : p = N-1 mod 4} (level theta = 1/2, Bombieri-Vinogradov):
  class-3 primes  < y  = N^(1/2 - delta)   (semi-linear),  level N^theta3, s3 = theta3/(1/2-delta)
  class-1 primes  < z0 = N^alpha           (semi-linear),  level N^theta1, s1 = theta1/alpha,  theta3+theta1 = 1/2
  (odd primes < z0 are all covered: class-1 by the second sieve, class-3 by the first)
Vector sieve:  S >= X V3(y) V1(z0) * L,  L = f3 F1 + F3 f1 - F3 F1.
Survivors that are NOT sums of two squares: n = q1 q2 m, q1<q2 = 3 mod 4, q1 >= y, m <= N^(2 delta) made of
class-1 primes >= z0.  Switching: sieve B = {N - q1 q2 m} for primes, linear upper sieve at level N^(1/2)
(BV for the convolution q1*q2*m):  #E <= 4 * Sing(N) * |B| / log N.
Normalising both sides by Sing(N) N / log^2 N:
   main  = e^-gamma L / (2 sqrt(alpha (1/2 - delta)))
   minus = 4 c,   c = |B| / (N / log N) = sum_m (1/m) * (1/4) * int_{1/2-delta}^{(1-mu)/2} du / (u (1 - mu - u)),  m = N^mu
"""
import numpy as np
exec(open(__file__.replace("net_yield.py", "sieve_funcs.py"), encoding="utf-8-sig").read().split("# ----")[0])
Fh, fh = SF[0.5]
EG = np.exp(np.euler_gamma)


def c_switch(delta, alpha, h=5e-4):
    # density of z0-rough class-1 m on the log scale mu in [0, 2 delta]: exp-convolution of (1/2) du/u on [alpha, 2delta]
    mu = np.arange(0, 2 * delta + h / 2, h)
    nu1 = np.where(mu >= alpha, 0.5 / np.maximum(mu, 1e-12), 0.0) * h
    total = np.zeros_like(mu); total[0] = 1.0           # m = 1
    term = total.copy()
    for j in range(1, int(2 * delta / alpha) + 1):       # m with j prime factors: nu1^{*j} / j!
        term = np.convolve(term, nu1)[: len(mu)] / j
        total += term
    lo = 0.5 - delta; hi = (1 - mu) / 2
    I = np.where(hi > lo, np.log(hi / lo) / (1 - mu) - np.log((1 - mu - hi) / (1 - mu - lo)) / (1 - mu), 0.0)
    # int du/(u(a-u)) = (1/a) [ln u - ln(a-u)]
    return float(np.sum(total * 0.25 * I))


def best(delta, alpha):
    t3 = np.linspace(0, 0.5, 501)
    s3, s1 = t3 / (0.5 - delta), (0.5 - t3) / alpha
    ok = (s3 > 0.05) & (s1 > 0.05) & (s3 < S_MAX) & (s1 < S_MAX)
    a, A_, b, B_ = at(fh, s3), at(Fh, s3), at(fh, s1), at(Fh, s1)
    L = np.where(ok, a * B_ + A_ * b - A_ * B_, -np.inf)
    i = int(np.argmax(L))
    main = np.exp(-np.euler_gamma) * L[i] / (2 * np.sqrt(alpha * (0.5 - delta)))
    minus = 4 * c_switch(delta, alpha)
    return main, minus, s3[i], s1[i]


print("sanity: c(delta, alpha) with m = 1 only (alpha > 2 delta) should be ~ delta for small delta:",
      c_switch(0.02, 0.05), "vs", 0.02)

results = []
for delta in np.arange(0.01, 0.245, 0.01):
    for alpha in np.arange(0.002, 0.2, 0.001):
        main, minus, s3, s1 = best(delta, alpha)
        if main - minus > 1e-3:
            results.append((int(np.ceil(1 / alpha)) - 1, delta, alpha, main, minus, s3, s1))

if not results:
    print("no (delta, alpha) gives a positive net yield")
else:
    results.sort(key=lambda r: (r[0], -(r[3] - r[4])))
    print(" k   delta  alpha    main    minus    net     s3     s1")
    seen = set()
    for r in results:
        if r[0] in seen or len(seen) >= 8:
            continue
        seen.add(r[0])
        k, d, a, m, mi, s3, s1 = r
        print(f"{k:3d}  {d:.3f}  {a:.3f}  {m:.4f}  {mi:.4f}  {m - mi:+.4f}  {s3:.3f}  {s1:.3f}")


