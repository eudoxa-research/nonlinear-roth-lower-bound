"""Exact certificate that the computed policy is optimal at m = 65 (or 85).
f(F) := exact value of the stationary policy pi from state F, as a rational combination of logs.
Check f(F) >= log T(F,F') + f(F')/2 for all F,F' by exact big-integer comparison.
Then W = f exactly (f >= B'f gives W <= f, and f = value of a feasible policy gives W >= f)."""
import sys, math, time, pickle
from fractions import Fraction
from collections import defaultdict
import numpy as np
m = int(sys.argv[1]) if len(sys.argv) > 1 else 65
Tfile = sys.argv[2] if len(sys.argv) > 2 else "T65.bin"
SQ = sorted({(x*x) % m for x in range(m)} - {0})
pairs = []; seen = set()
for s in SQ:
    if s in seen: continue
    pairs.append((s, (-s) % m)); seen |= {s, (-s) % m}
NP = len(pairs); NS = 1 << NP; S = NS - 1
T = np.fromfile(Tfile, dtype=np.uint8).reshape(NS, NS).astype(np.int64)
logT = np.log(T.astype(float))
# float value iteration to get a candidate policy
W = np.zeros(NS)
for it in range(500):
    Wn = (logT + 0.5*W[None, :]).max(axis=1)
    if np.abs(Wn - W).max() < 1e-14: W = Wn; break
    W = Wn
pol = (logT + 0.5*W[None, :]).argmax(axis=1)
print("float W(S) =", repr(W[S]))
# exact policy values: f(F) = sum_n c_n log n
def policy_value(F):
    path = []; idx = {}
    while F not in idx:
        idx[F] = len(path); path.append(F); F = pol[F]
    pre = idx[F]; per = len(path) - pre
    c = defaultdict(Fraction)
    for i in range(pre):
        c[int(T[path[i], pol[path[i]]])] += Fraction(1, 2**i)
    g = Fraction(1, 1) / (1 - Fraction(1, 2**per))
    for j in range(per):
        i = pre + j
        c[int(T[path[i], pol[path[i]]])] += g / 2**i
    return dict(c), pre, per
f = [None]*NS; pres = []; pers = []
t0 = time.time()
for F in range(NS):
    f[F], pre, per = policy_value(F)
    pres.append(pre); pers.append(per)
print("policy values computed in", time.time()-t0, "s; max pre-period", max(pres), "max period", max(pers), "periods:", sorted(set(pers)))
def fval(c): return sum(float(v)*math.log(n) for n, v in c.items())
print("exact f(S) =", f[S], " numeric", fval(f[S]))
# exact comparison: f(F) - f(F')/2 - log T >= 0  <=>  prod n^(D e_n) >= T^D  with D common denominator
def compare(cF, cFp, t):
    e = defaultdict(Fraction)
    for n, v in cF.items(): e[n] += v
    for n, v in cFp.items(): e[n] -= v/2
    e[t] -= 1
    D = 1
    for v in e.values(): D = D*v.denominator//math.gcd(D, v.denominator)
    lhs = 1; rhs = 1
    for n, v in e.items():
        k = int(v*D)
        if k > 0: lhs *= n**k
        elif k < 0: rhs *= n**(-k)
    return (lhs > rhs) - (lhs < rhs)   # 1: strict, 0: equality, -1: violated
t0 = time.time(); viol = 0; eq = defaultdict(list); nbig = 0
fv = np.array([fval(c) for c in f])
Dm = fv[:, None] - 0.5*fv[None, :] - logT      # float screen, error < 1e-12
cand = np.argwhere(Dm <= 1e-6)
print("float minimum of f(F) - f(F')/2 - log T:", Dm.min(), " pairs needing exact check:", len(cand))
for F, Fp in cand:
    F = int(F); Fp = int(Fp); t = int(T[F, Fp]); nbig += 1
    r = compare(f[F], f[Fp], t)
    if r < 0:
        viol += 1; print("VIOLATION", F, Fp, t, Dm[F, Fp])
    elif r == 0:
        eq[F].append((Fp, t))
print("exact checks needed:", nbig, " violations:", viol, " time", time.time()-t0)
# optimal actions (exact equalities) reachable from S
seen = set(); frontier = [S]; sizes_seen = set()
while frontier:
    F = frontier.pop()
    if F in seen: continue
    seen.add(F)
    acts = eq[F]
    print("state pairs", bin(F).count('1'), "exact optimal actions (|R|, next pairs):", sorted({(t, bin(Fp).count('1')) for Fp, t in acts}))
    for Fp, t in acts: frontier.append(Fp)
print("states on optimal paths:", len(seen))
pickle.dump({"f": f, "eq": dict(eq), "pol": pol}, open(f"cert{m}.pkl", "wb"))
