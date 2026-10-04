"""
Rigorous (outward-rounded) check of the MAIN-TERM inequality for the fixed parameters
    delta = 0.10, alpha = 0.007, lam = 0.15, beta1 = 0.30, K = 4
    net = (L - lam*R) / norm - K*c  > 0,   norm = 2 e^gamma sqrt(alpha (1/2 - delta))
This certifies the numerical inequality only.  It does NOT certify the sieve lemmas (L1-L6) or the o(1) terms.

Assumptions used (standard facts about the kappa = 1/2, beta = 1 sieve functions):
  F decreasing, f increasing on s > 0;  F(s) = A/sqrt(s) on (0,2];  f(s) = A arccosh(sqrt s)/sqrt s on [1,3];
  (sqrt(s) F)' = f(s-1) / (2 sqrt s),  (sqrt(s) f)' = F(s-1) / (2 sqrt s).
Rounding: closed forms and logs via mpmath interval arithmetic; the long recurrence in IEEE doubles with one
nextafter outward step after every operation (+,-,*,/,sqrt are correctly rounded in IEEE 754).
"""
import math
import numpy as np
from mpmath import iv, mp
iv.prec = 80
mp.prec = 80

dn = lambda x: math.nextafter(x, -math.inf)
up = lambda x: math.nextafter(x, math.inf)

A_iv = 2 * iv.sqrt(iv.exp(iv.euler) / iv.pi)
A_lo, A_hi = float(A_iv.a), float(A_iv.b)
A_lo, A_hi = dn(A_lo), up(A_hi)

H_DEN = 1000                      # grid step h = 1/1000, grid point i <-> s = i/1000 (exact rationals)
S_TOP = 80
n = S_TOP * H_DEN
Flo = np.zeros(n + 1); Fhi = np.zeros(n + 1); flo = np.zeros(n + 1); fhi = np.zeros(n + 1)
Flo[0] = Fhi[0] = math.inf
for i in range(1, 3 * H_DEN + 1):
    s = iv.mpf([i, i]) / H_DEN
    if i <= 2 * H_DEN:
        v = A_iv / iv.sqrt(s); Flo[i], Fhi[i] = dn(float(v.a)), up(float(v.b))
    if i >= H_DEN:
        v = A_iv * iv.log(iv.sqrt(s) + iv.sqrt(s - 1)) / iv.sqrt(s); flo[i], fhi[i] = max(0.0, dn(float(v.a))), up(float(v.b))

sq_lo = [0.0] * (n + 1); sq_hi = [0.0] * (n + 1)          # enclosures of sqrt(i/1000)
for i in range(n + 1):
    v = iv.sqrt(iv.mpf([i, i]) / H_DEN); sq_lo[i], sq_hi[i] = dn(float(v.a)), up(float(v.b))

gF_lo, gF_hi = dn(sq_lo[2 * H_DEN] * Flo[2 * H_DEN]), up(sq_hi[2 * H_DEN] * Fhi[2 * H_DEN])
gf_lo, gf_hi = dn(sq_lo[3 * H_DEN] * flo[3 * H_DEN]), up(sq_hi[3 * H_DEN] * fhi[3 * H_DEN])
for i in range(2 * H_DEN, n):
    j = i - H_DEN                                            # lagged index: t-1 in [j, j+1] (grid units)
    w_lo, w_hi = dn(sq_lo[i + 1] - sq_hi[i]), up(sq_hi[i + 1] - sq_lo[i])   # int_{s_i}^{s_i+1} dt/(2 sqrt t)
    # F step: f(t-1) in [flo[j], fhi[j+1]] since f increasing
    gF_lo = dn(gF_lo + dn(flo[j] * w_lo)); gF_hi = up(gF_hi + up(fhi[j + 1] * w_hi))
    Flo[i + 1] = max(Flo[i + 1], dn(gF_lo / sq_hi[i + 1])) if i + 1 <= 2 * H_DEN else dn(gF_lo / sq_hi[i + 1])
    Fhi[i + 1] = up(gF_hi / sq_lo[i + 1])
    if i >= 3 * H_DEN:
        # f step: F(t-1) in [Flo[j+1], Fhi[j]] since F decreasing
        gf_lo = dn(gf_lo + dn(Flo[j + 1] * w_lo)); gf_hi = up(gf_hi + up(Fhi[j] * w_hi))
        flo[i + 1] = dn(gf_lo / sq_hi[i + 1]); fhi[i + 1] = up(gf_hi / sq_lo[i + 1])


def F_up(s):   # upper bound of F at real s > 0 (F decreasing -> use grid point at or below s)
    if s <= 2:
        return up(A_hi / dn(math.sqrt(s)))
    i = min(int(math.floor(s * H_DEN)) - 1, n)
    return Fhi[i] if s * H_DEN <= n + 1 else Fhi[n]


def F_dn(s):   # lower bound of F (use grid point at or above s; F >= 1 beyond the grid is NOT assumed)
    i = int(math.ceil(s * H_DEN)) + 1
    return Flo[i] if i <= n else 0.0


def f_dn(s):   # lower bound of f (f increasing -> grid point at or below s)
    i = int(math.floor(s * H_DEN)) - 1
    return flo[min(i, n)] if i >= 0 else 0.0


print(f"sieve functions: F(80) in [{Flo[n]:.6f}, {Fhi[n]:.6f}],  f(80) in [{flo[n]:.6f}, {fhi[n]:.6f}]")
print(f"                 F(2.5) in [{Flo[2500]:.6f}, {Fhi[2500]:.6f}],  f(2.5) in [{flo[2500]:.6f}, {fhi[2500]:.6f}]")

delta, alpha, lam, beta1, K = 0.10, 0.007, 0.15, 0.30, 4.0      # all exactly representable as decimals below
D, AL, B1 = mp.mpf('0.10'), mp.mpf('0.007'), mp.mpf('0.30')
# exact rationals: 1/lam + 1/beta1 = 20/3 + 10/3 = 10, so w_n > 0 forces Omega(n) < 10, i.e. k = 9

# ---- L: lower bound at a fixed split theta3 (any split is admissible) -------------------------------
best = None
for t3 in np.arange(0.40, 0.4999, 0.0005):
    s3, s1 = t3 / (0.5 - delta), (0.5 - t3) / alpha
    Lv = f_dn(s3) * F_dn(s1) + F_up(s3) * (f_dn(s1) - F_up(s1))
    if best is None or Lv > best[0]:
        best = (Lv, t3)
t3 = best[1]
s3, s1 = t3 / (0.5 - delta), (0.5 - t3) / alpha
# L = f3 F1 + F3 (f1 - F1), with f1 - F1 <= 0:  L >= f3_lo F1_lo + F3_hi (f1_lo - F1_hi)
L_lo = dn(dn(f_dn(s3) * F_dn(s1)) + dn(F_up(s3) * dn(f_dn(s1) - F_up(s1))))

# ---- R: upper bound,  R = 1/2 int_alpha^beta1 (1 - u/beta1) U(u) du/u,  U <= F(s3') F(s1'(u)) -----------
NU = 6000
R_hi = 0.0
for k in range(NU):
    ua = alpha + (beta1 - alpha) * k / NU; ub = alpha + (beta1 - alpha) * (k + 1) / NU
    ua, ub = dn(ua), up(ub)
    lev = 0.5 - ub                                    # remaining level at the worst end of the cell
    bestU = math.inf
    for frac in np.linspace(0.70, 0.995, 60):         # split chosen per cell; any split gives a valid upper bound
        t = lev * frac
        # s1 = (0.5 - u - t)/alpha is smallest at u = ub -> F(s1) largest there
        val = F_up(t / (0.5 - delta)) * F_up(dn((0.5 - ub - t) / alpha))
        bestU = min(bestU, val)
    weight = up(up(1 - ua / beta1) / ua)              # (1 - u/beta1)/u is decreasing: max at ua
    R_hi = up(R_hi + up(up(weight * bestU) * up(ub - ua)))
R_hi = up(0.5 * R_hi)

# ---- c: upper bound (masses of class-1 m shifted LEFT, I(mu) decreasing) -----------------------------------
# I(mu) = ln((1/2+delta-mu)/(1/2-delta)) / (1-mu) = ln((0.6-mu)/0.4)/(1-mu); its derivative has the sign of
# -1/(0.6-mu) + I(mu) <= -1/0.6 + ln(1.5)/0.8 < 0, so I is decreasing on [0, 0.2] and left-shifting is an upper bound.
# Cell masses are the exact integrals (1/2) ln((k+1)/k) (the float scan used the cruder 0.5 h / mu, hence its larger c).
HM = mp.mpf(1) / 2000                                  # mu-grid 0.0005: alpha = 14 HM, 2 delta = 400 HM
i_al, i_top = 14, 400
mass = np.zeros(i_top + 1)
for k in range(i_al, i_top):
    v = iv.mpf(0.5) * iv.log(iv.mpf([k + 1, k + 1]) / k)    # (1/2) int_{k h}^{(k+1) h} du/u, placed at k h
    mass[k] = up(float(v.b))
total = np.zeros(i_top + 1); total[0] = 1.0; term = total.copy()
for j in range(1, i_top // i_al + 1):
    term = np.convolve(term, mass)[: i_top + 1] / j
    total += term
total *= (1 + 1e-9)                                    # covers float rounding in the convolution (<1e-12 rel.)
c_hi = 0.0
for k in range(i_top + 1):
    if total[k] == 0:
        continue
    mu = iv.mpf([k, k]) * HM
    lo_, hi_ = iv.mpf('0.5') - iv.mpf('0.1'), (1 - mu) / 2
    if hi_.a <= lo_.b:
        continue
    a = 1 - mu
    I = (iv.log(hi_ / lo_) - iv.log((a - hi_) / (a - lo_))) / a
    c_hi = up(c_hi + up(total[k] * up(0.25 * float(I.b))))

# ---- assemble ---------------------------------------------------------------------------------------------
norm = 2 * iv.exp(iv.euler) * iv.sqrt(iv.mpf('0.007') * (iv.mpf('0.5') - iv.mpf('0.1')))
norm_lo, norm_hi = dn(float(norm.a)), up(float(norm.b))
num_lo = dn(L_lo - up(lam * R_hi))
net_lo = dn((dn(num_lo / norm_hi) if num_lo >= 0 else dn(num_lo / norm_lo)) - up(K * c_hi))
print(f"theta3 = {t3:.4f}  (s3 = {s3:.4f}, s1 = {s1:.3f})")
print(f"L  >= {L_lo:.6f}")
print(f"R  <= {R_hi:.6f}")
print(f"c  <= {c_hi:.6f}")
print(f"norm in [{norm_lo:.6f}, {norm_hi:.6f}]")
print(f"certified:  net >= {net_lo:+.6f}   ->  {'POSITIVE' if net_lo > 0 else 'NOT certified'}")
