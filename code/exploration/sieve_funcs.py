"""
Beta-sieve functions f_k, F_k for kappa = 1/2 (beta = 1, Iwaniec semi-linear) and kappa = 1 (beta = 2, Rosser-Iwaniec linear).
  (s^k F(s))' = k s^(k-1) f(s-1),  (s^k f(s))' = k s^(k-1) F(s-1)
  F(s) = A s^-k  for 0 < s <= beta+1,   f(s) = 0 for s <= beta.
The constant A is not hard-coded: it is fitted so that F(s), f(s) -> 1 as s -> infinity,
then compared with the known closed forms (A = 2e^gamma for kappa=1, A = 2 sqrt(e^gamma/pi) for kappa=1/2).
"""
import numpy as np

H = 1e-4
S_MAX = 30.0
grid = np.arange(0, S_MAX + H / 2, H)
n1 = int(round(1 / H))


def solve(kappa, beta, A=1.0):
    F = np.zeros_like(grid); f = np.zeros_like(grid)
    i0F = int(round((beta + 1) / H)); i0f = int(round(beta / H))
    F[1:i0F + 1] = A * grid[1:i0F + 1] ** -kappa
    gF = grid ** kappa * F; gf = np.zeros_like(grid)  # g = s^kappa * (function)
    if kappa == 0.5:
        # on 1 <= s <= 3, F(s-1) = A (s-1)^-1/2 is singular at s=1: use the exact integral
        # sqrt(s) f(s) = (A/2) int_1^s dt / sqrt(t(t-1)) = A arccosh(sqrt s)
        j = int(round(3 / H))
        gf[i0f:j + 1] = A * np.arccosh(np.sqrt(grid[i0f:j + 1]))
        f[i0f:j + 1] = gf[i0f:j + 1] / grid[i0f:j + 1] ** kappa
        i0f = j
    # trapezoid integration; f(s-1), F(s-1) are always already known (lag 1 >= step)
    for i in range(min(i0f, i0F) + 1, len(grid)):
        s0, s1 = grid[i - 1], grid[i]
        if i > i0F:
            d0 = kappa * s0 ** (kappa - 1) * f[i - 1 - n1]; d1 = kappa * s1 ** (kappa - 1) * f[i - n1]
            gF[i] = gF[i - 1] + H * (d0 + d1) / 2; F[i] = gF[i] / s1 ** kappa
        if i > i0f:
            d0 = kappa * s0 ** (kappa - 1) * F[i - 1 - n1]; d1 = kappa * s1 ** (kappa - 1) * F[i - n1]
            gf[i] = gf[i - 1] + H * (d0 + d1) / 2; f[i] = gf[i] / s1 ** kappa
    return F, f


def at(arr, s):
    return np.interp(s, grid, arr)


EG = np.exp(np.euler_gamma)
SF = {}
for kappa, beta, A_known in ((1.0, 2, 2 * EG), (0.5, 1, 2 * np.sqrt(EG / np.pi))):
    F1, f1 = solve(kappa, beta, 1.0)          # linear in A, so solve with A=1 and rescale
    A_fit = 2 / (at(F1, S_MAX) + at(f1, S_MAX))
    F, f = F1 * A_known, f1 * A_known
    SF[kappa] = (F, f)
    print(f"kappa={kappa}: A fitted from limit = {A_fit:.5f}, known closed form = {A_known:.5f}; "
          f"F(30)={at(F,30):.6f}, f(30)={at(f,30):.6f}")
    for s in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0):
        print(f"   s={s:3.1f}: F={at(F,s):.4f}  f={at(f,s):.4f}  f/F={at(f,s)/at(F,s):.4f}")

# closed-form checks
F, f = SF[1.0]
print("linear check f(3) vs 2e^g log2/3:", at(f, 3), 2 * EG * np.log(2) / 3)
F, f = SF[0.5]
print("semi-linear check f(2) vs sqrt(e^g/(2pi))*2*arccosh(sqrt2):",
      at(f, 2), np.sqrt(EG / (np.pi * 2)) * 2 * np.arccosh(np.sqrt(2)))

# ---------------------------------------------------------------------------
# Vector sieve for n = N - p, p = N-1 mod 4 (so n = 1 mod 4):
#   class-3 primes sifted up to y = N^eta (semi-linear),  class-1 primes up to z = N^(1/v) (semi-linear).
#   level split theta3 + theta1 = theta; s3 = theta3/eta, s1 = theta1*v.
#   Brudern-Fouvry: ab >= a- b+ + a+ b- - a+ b+  ->  L = f3 F1 + F3 f1 - F3 F1.
#   Survivors: no class-1 factor < z, no class-3 factor < y.  If eta >= 1/2 (and n = 1 mod 4)
#   there is no class-3 factor at all, so n is a sum of two squares with Omega(n) < v.
# ---------------------------------------------------------------------------
Fh, fh = SF[0.5]


def best_L(theta, eta, v, steps=400):
    best = -np.inf
    for t3 in np.linspace(0, theta, steps):
        s3, s1 = t3 / eta, (theta - t3) * v
        if s3 > S_MAX or s1 > S_MAX:
            continue
        a, A_ = at(fh, s3), at(Fh, s3); b, B_ = at(fh, s1), at(Fh, s1)
        best = max(best, a * B_ + A_ * b - A_ * B_)
    return best


print("\n(A) BV level theta = 1/2, eta = 1/2 (no class-3 factor can survive):")
for v in (5, 10, 20, 50, 100):
    print(f"   v={v:3d}: max L = {best_L(0.5, 0.5, v):+.4f}")

print("\n(B) hypothetical level theta > 1/2 (e.g. BFI-type), eta = 1/2: smallest k = ceil(v)-1 with L > 0")
for theta in (0.55, 4 / 7, 0.6, 2 / 3, 0.75, 1.0):
    ks = next(([v] for v in np.arange(2.0, 200.0, 0.25) if best_L(theta, 0.5, v, 200) > 1e-3), [])
    if ks:
        v = ks[0]; print(f"   theta={theta:.4f}: L>0 first at v={v:.2f}  ->  k = {int(np.ceil(v)) - 1}")
    else:
        print(f"   theta={theta:.4f}: no v < 200 works")

