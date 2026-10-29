Ancillary files for "An improved lower bound in the nonlinear Roth theorem" (revised 13 September 2026).

verify.py      Appendix A. Checks the chain at m = 145 behind Theorem 1.1 (sizes 10, 28 then 20, 21, 20, 22 periodic)
               and the chain (7, 17, 10, 15) at m = 65 behind Theorem 1.2, and recomputes both exponents from the sizes
               through Corollary 2.3.  Standard library only.                      python3 verify.py
verify85.py    Checks an explicit optimal chain at m = 85 (three sets and a four-cycle) and recomputes Gamma(85).
chain145.json  The six sets of Proposition 3.1 as JSON.

dp_general.c   Section 4. Builds the table T(F,F') of clique numbers (Lemma 4.2) for a squarefree modulus m <= 127 with
               -1 a square mod m, runs value iteration, prints W(S), Gamma(m) and the optimal policy, writes the table.
                   gcc -O3 -o dp_general dp_general.c -lm ; ./dp_general 65 T65.bin (2 s) ; ./dp_general 85 T85.bin (4 min)
exact_cert.py  Section 4.2. Exact big-integer certificate that the computed policy is optimal (needs numpy).
                   python3 exact_cert.py 65 T65.bin ; python3 exact_cert.py 85 T85.bin

cliq.c         Section 5. Maximum cliques in Cayley graphs on Z/mZ with a symmetric forbidden difference set, m <= 1024,
               multi-word bitsets, greedy-colouring bound, vertex 0 fixed; modes: omega only, omega plus one clique,
               enumeration of up to K maximum cliques containing 0; optional node budget (truncated searches still
               return valid cliques, so the restricted program stays a lower bound).
                   gcc -O3 -o cliq cliq.c
rdp.py         Section 5. The restricted (lower-bound) dynamic program.  Needs ./cliq in the working directory.
                   python3 rdp.py m rounds enumcap topn seed minpairs ktop budget twopair
               e.g.  python3 rdp.py 145 16 5000 20 1 8 80 0 0   reaches 0.7750918639 in a few minutes and writes best_145.pkl
               (the chain of Proposition 3.1 came from such a run; different seeds and parameters explore different states).
