import math
# Explicit optimal chain at m = 85 (Proposition 4.3): R1, R2, R3, then the four-cycle R4, R5, R6, R7, R4, ...
m = 85
R = [
    [37, 43, 45, 74, 76, 82, 84],
    [17, 19, 22, 24, 42, 44, 47, 49, 52, 54, 72, 74, 77, 79, 82, 84],
    [10, 13, 16, 27, 30, 33, 44, 47, 50, 61, 64, 67, 78, 81, 84],
    [28, 38, 39, 40, 41, 42, 52, 70, 80, 81, 82, 83, 84],
    [10, 13, 16, 24, 27, 30, 41, 44, 47, 54, 61, 64, 71, 78, 81, 84],
    [60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84],
    [10, 13, 16, 27, 30, 33, 44, 47, 50, 61, 64, 67, 78, 81, 84]
]
diffs = lambda X: {(a - b) % m for a in X for b in X if a != b}
sqim = lambda X: {((a - b) ** 2) % m for a in X for b in X if a != b}
SQ = {(x * x) % m for x in range(m)} - {0}
for X in R:
    assert len(set(X)) == len(X) and all(0 <= x < m for x in X)
assert diffs(R[0]) & SQ == set()                       # condition at n = 0
for k in range(6):
    assert diffs(R[k + 1]) & sqim(R[k]) == set()       # chain conditions
assert diffs(R[3]) & sqim(R[6]) == set()               # R8 = R4 closes the cycle
sizes = [len(X) for X in R]                            # [7, 16, 15, 13, 16, 13, 15]
pre, per = 3, 4
W = sum(math.log(sizes[i]) / 2 ** i for i in range(pre))
W += sum(math.log(sizes[pre + j]) / 2 ** (pre + j) for j in range(per)) / (1 - 2.0 ** -per)
gamma = 0.5 + W / (4 * math.log(m))
closed = 0.5 + (30*math.log(7) + 17*math.log(16) + 8*math.log(15) + 5*math.log(13)) / (120*math.log(85))
print("verified; sizes =", sizes, "; Gamma(85) =", gamma, "; closed form =", closed)
