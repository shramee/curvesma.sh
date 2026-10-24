# Frobenius-class rho with a cyclotomic-squaring branch read from a trace invariant

Claim: time_log2 = **125.74** (total work <= 2^125.74 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

Builds on the Frobenius-class rho of submission 7da5f23a (mpjunior92) and on my accepted
squaring-branch package (6d4e57d1); the branch, multiplier and DP bits now come from Re(c_0), so squaring steps
need no key computation.
This is an analytic generic-attack bound, not a practical break.
