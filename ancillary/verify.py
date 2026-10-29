import math

def check(m, A, B):
    """A = pre-period sets, B = cycle sets (Corollary 2.3); returns gamma."""
    C = A + B
    for X in C:
        assert len(set(X)) == len(X) and all(0 <= x < m for x in X)
    SQ = {(x * x) % m for x in range(m)} - {0}
    diffs = lambda X: {(a - b) % m for a in X for b in X if a != b}
    sqim = lambda X: {((a - b) ** 2) % m for a in X for b in X if a != b}
    assert diffs(C[0]) & SQ == set()                   # condition (1)
    for X, Y in zip(C, C[1:]):
        assert diffs(Y) & sqim(X) == set()             # condition (2)
    assert diffs(B[0]) & sqim(B[-1]) == set()          # wrap-around
    q, p = len(A), len(B)
    g = 0.5 + sum(math.log(len(X), m) / 2 ** (i + 2) for i, X in enumerate(A))
    cyc = sum(math.log(len(X), m) / 2 ** (q + j + 2) for j, X in enumerate(B))
    g += cyc / (1 - 2.0 ** -p)
    print("m =", m, [len(X) for X in A], [len(X) for X in B], "gamma =", g)
    return g

A1 = [0, 48, 55, 62, 69, 76, 117, 124, 131, 138]
A2 = [0, 10, 12, 22, 25, 27, 35, 37, 50, 52, 60, 62, 65, 67, 75, 77,
      90, 92, 100, 102, 105, 115, 117, 127, 130, 132, 140, 142]
B1 = [0, 15, 18, 26, 29, 44, 47, 55, 58, 73, 76, 84, 87, 102, 105, 113,
      116, 131, 134, 142]
B2 = [0, 7, 15, 38, 46, 54, 61, 69, 76, 77, 84, 91, 92, 99, 107, 114,
      122, 129, 130, 137, 144]
B3 = [0, 9, 13, 22, 31, 35, 44, 53, 57, 66, 75, 79, 88, 97, 101, 110,
      114, 123, 132, 136]
B4 = [0, 22, 25, 32, 35, 47, 50, 57, 60, 72, 75, 82, 85, 92, 95, 107,
      110, 117, 120, 132, 135, 142]
g145 = check(145, [A1, A2], [B1, B2, B3, B4])           # 0.7750918639366062

C1 = [5, 11, 16, 22, 53, 59, 64]
C2 = [0, 1, 6, 11, 16, 17, 22, 27, 32, 33, 38, 43, 48, 49, 54, 59, 64]
C3 = [1, 6, 11, 44, 49, 54, 56, 59, 61, 64]
C4 = [1, 6, 12, 14, 19, 25, 27, 32, 38, 40, 45, 51, 53, 58, 64]
g65 = check(65, [], [C1, C2, C3, C4])                   # 0.7732006139958438
