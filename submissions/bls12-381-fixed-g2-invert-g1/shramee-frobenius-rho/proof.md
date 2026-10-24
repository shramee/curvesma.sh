# Frobenius-class parallel rho with a cyclotomic-squaring branch read from a trace invariant

**Claim.** For uniform nonidentity Q in G2 and uniform z in GT the algorithm below
recovers P in G1 with e(P,Q) = z with probability at least 0.39, using at most
**2^125.74** dense-Fp12-multiplication equivalents (C = 20,000), including setup, the two
pairings, table construction, walk starts, distinguished-point (DP) handling, cycle
escapes, the DP detection delay, output and verification. Peak memory is below 2^73
bytes, preprocessing (included) at most 2^31 units, nonuniform advice 0 bytes. The success
bound rests on the usual random-mapping heuristic for rho walks, adjusted for the branch
mix (section 4), and is supported by toy-field experiments and by running the production
step rule in the real GT. No production-size run is claimed.

The class construction (Frobenius classes of GT, the 4-multiplication key, the twisted
multiplier table that never applies a Frobenius to the walk element) is the one published
in submission 7da5f23a-6633-46b1-9c03-11f6cfedeeec (mpjunior92), which this package builds
on. New here: (i) a cyclotomic-squaring branch that is cheaper than the table multiply;
(ii) a fresh instruction ledger; (iii) an analysis of how the branch mix changes the
collision law and the fruitless cycles; (iv) execution of the production step rule in the
real BLS12-381 GT; (v) explicit memory and report caps; (vi) the branch, multiplier and DP bits are read from the free class
invariant T = Re(c_0) instead of the minimum key K, so a squaring step needs no key computation at all.
The package `frobenius-rho-squaring-mix` (my earlier submission 6d4e57d1, accepted as valid at 125.80)
is the same design with the bits read from K; this one differs only by (vi) and its ledger.

## 1. Reduction

Compute g = e(P0, Q) with the pinned convention pairing(Q, P0). Q has prime order r and
the pairing is nondegenerate, so g != 1 generates GT (prime order r). z = g^x for a unique
x in [0, r). Find x, output P = [x]P0, check e(P, Q) = z. If z = 1 output infinity. The
answer is for the supplied Q.

## 2. Arithmetic facts (checked in `experiment.py`, parts A and C)

Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u (the blst tower). With
u = w^6 - 1 this is the profile field, w^12 - 2w^6 + 2 = 0, and conversion of coefficients
is free: c_m = a_m - b_m, c_{m+6} = b_m.

**F1 Frobenius classes.** (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}.
p has order exactly 12 modulo r, p^6 = -1 and 1 + p^4 + p^8 = 0 (mod r), so the maps
y -> y^{p^k} form a group of order 12 acting on GT, containing inversion. Every y != 1
lies in a class [y] of exactly 12 elements; N = (r-1)/12 classes, log2 N = 251.2721.

**F2 Keys.** Put zeta = xi^{(p^2-1)/6}; it lies in Fp and is a primitive 6th root of unity
(zeta^2 = zeta - 1). Let c_1 = a + b u be the w^1 coefficient of y and B = Re(conj(c_1) d_1)
= a Re(d_1) + b Im(d_1). The w^1 coefficient of y^{p^(k+1)} is conj(c^{(k)}) d_1, and
conj(d_1) d_1 = zeta, hence
`key_k := Re(w^1 coeff of y^{p^k}) = a zeta^j (k = 2j), B zeta^j (k = 2j+1)`.
With zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta} all 12 keys cost **4 Fp multiplications**
(a zeta, a Re d_1, b Im d_1, B zeta) plus 8 subtractions/negations and one addition. The
12 keys are distinct unless two of them coincide by a linear accident (probability <= 66/p
per element; ties are charged as failure). K([y]) = min_k key_k is a class invariant,
kappa(y) = argmin_k key_k, and y^{p^kappa} is the same element for all 12 members of the
class (the canonical representative y^). K is the minimum of 12 values that behave like
independent uniform residues (they come in negation pairs, so K lies in [0, p/2]); its low 200
bits are uniform up to statistical distance about 2^200/p = 2^-181 (the density of the minimum is smooth).

**F2b Trace invariant.** The w^0 coefficient of y^{p^k} is conj^k(c_0) (no twist), so T([y]) := Re(c_0) is the same
for all 12 conjugates, and 12 Re(c_0) is the trace of y to Fp. T is read directly from the coefficients (no
multiplication) and, for distinct classes, is an essentially uniform residue (collisions of T between distinct classes
have probability about 1/p per pair). It plays the role of the class identity for the bits below; K is needed only to
find kappa. `experiment.py` part A checks that Re(c_0) is equal over all 12 conjugates of random GT elements.

**F3 Cyclotomic squaring (Granger-Scott).** GT lies in the cyclotomic subgroup, where
squaring is cheaper than a product. Write y = A + B w + C w^2 with A=(c0,c3), B=(c1,c4),
C=(c2,c5) in Fp4 = Fp2[v]/(v^2-xi), v = w^3, and (A0,A1) = A^2, (B0,B1) = B^2,
(C0,C1) = C^2 computed in Fp4 with 3 Fp2 squarings each. Then
`c0' = 3A0 - 2c0, c3' = 3A1 + 2c3, c1' = 3 xi C1 + 2c1, c4' = 3C0 - 2c4, c2' = 3B0 - 2c2, c5' = 3B1 + 2c5.`
`experiment.py` verifies y^2 = formula on real BLS12-381 GT elements and counts the
operations of the implementation: **9 Fp2 squarings and 34 Fp2 additions/subtractions**
(the multiplications by xi are counted as one Fp2 addition).

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 *
2529403 * 52437899 * 254760293^2 (verified by multiplication in `experiment.py`). The order
of 2 modulo r is (r-1)/16 (about 2^250.9), so 2^t lies in the order-12 group <p> only if t is a
multiple of about 2^247.3.

## 3. Algorithm

Parameters: R = 1024 table multipliers, squaring probability s = 3/8, theta = 2^-60 (DP
density), M = 2^20 parallel walks. Bits of the class invariant T = Re(c_0) (read as an integer in [0,p)): T mod 8 selects
the branch (squaring iff T mod 8 < 3), h = (T >> 3) mod 1024 selects the multiplier, the class
is distinguished iff (T >> 13) mod 2^60 = 0. The fingerprint of a class used for matching reports is
F = (T >> 73) mod 2^128, the 128 bits above the branch, multiplier and DP bits, so they are independent.

Precompute M_h = g^{a_h} z^{b_h} for h < R (fresh uniform a_h, b_h mod r) and
T[h][t] = M_h^{p^t}, t = 0..11, with the scalar tables a_h p^t, b_h p^t mod r
(12,288 entries of 576 bytes, 7.1 MB).

A walk keeps some member y of the current class together with (alpha, beta) mod r such that
y = g^alpha z^beta. Each walk starts at fresh uniform (alpha, beta). One step:

```
1  T <- Re(c_0) of y (no multiplication); branch, h and DP bits from T
2  if class distinguished: compute the 12 keys (4 Fp mult) -> kappa; report (F, p^kappa*alpha, p^kappa*beta)
   into a chained multimap keyed by F
3  2-cycle / checkpoint tests on T (section 3.3)
4  if T mod 8 < 3 (squaring branch):  y <- y^2 (Granger-Scott);  alpha,beta <- 2*alpha, 2*beta
   [no keys, no twist]
   else (multiply branch): compute the 12 keys (4 Fp mult) -> kappa; h <- (T>>3) mod 1024;
        t <- (12 - kappa) mod 12;  y <- y * T[h][t]   (ONE dense Fp12 multiplication);
        alpha += a_h p^t;  beta += b_h p^t
```

**3.1 The step is a map on classes.** Write y^ = y^{p^kappa}. Multiply branch:
y * M_h^{p^{-kappa}} = (y^ * M_h)^{p^{-kappa}} lies in the class [y^ M_h]. Squaring branch:
y^2 = (y^^2)^{p^{-kappa}} lies in [y^^2], because squaring commutes with Frobenius. Both
depend only on the class (T is a class function; kappa and y^ enter only through the conjugation
absorbed by the twist), so walks that reach the same class stay together and reach the same DP.

**3.2 Collision.** The multimap keeps every report with a given F (no entry is ever replaced), so a
useful class cannot be shadowed. Two reports with equal F from different walks are the same class, up to
a false match of probability 2^-128 per pair of distinct classes: with D <= 2^65.95 reports there are about
D^2/2^129 = 2^2.7 false candidates in expectation; every candidate is verified (one modular inverse and
g^x = z, <= 600 units), which is negligible. Their canonical exponents (alpha^, beta^) = p^kappa (alpha, beta) satisfy
alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1).
beta^1 = beta^2 happens with probability 1/r. Output P = [x]P0 after checking g^x = z.

**3.3 Fruitless cycles.** A closed walk of the class graph is *useless* if the canonical
exponent pair returns unchanged. After l steps the pair is Lambda (alpha^_0, beta^_0) + tau
where Lambda = 2^t p^m (t squaring steps, m a sum of twists) and tau = sum_h c_h (a_h, b_h)
with c_h in Z[p]. Because (alpha^_0, beta^_0) is uniform and independent of the table, the pair
returns unchanged with probability about 1/r unless Lambda = 1 and every c_h = 0.
* **No squaring step can occur in a useless cycle.** Lambda = 1 with t >= 1 needs
  2^t in <p> (F4), i.e. t >= 2^247.
* A useless cycle made of multiply steps needs every table index used to occur at least twice
  with coefficients p^{e_1} + p^{e_2} + ... = 0 in Z[p]/(Phi_12(p)): for two occurrences
  e_2 = e_1 + 6 (probability 1/12), for three {e, e+4, e+8}, and so on; for any number of occurrences
  the probability is at most 1/12 (the last twist must hit one value). For 2-cycles this gives the
  model rate (5/8)^2/(12 R) <= 1/(12 R) = 2^-13.6 per step (measured 4.0e-4 to 5.6e-4 per step across
  seeds versus 5.1e-4 predicted at R = 64, section 7). For l = 2k indices drawn independently and
  uniformly (heuristic H2) the probability is at most (2k)!/(2^k k!) (12R)^-k (pairings dominate), which
  is below 2^-288 for 65 <= l <= 147. For l >= 148 all l consecutive steps must be multiply steps, which
  has probability (5/8)^l < 2^-100 (the branch bits are bits of T, H2). Such cycles are never met in
  practice; if one is, the stuck-walk rule restarts the walk after 20/theta = 2^64.3 steps, which adds
  less than 2^-30 relative work at that rate.
* **Detection.** The step before last is remembered; if its T equals the current T the two
  classes form a 2-cycle and the walk replaces the representative of the class with the smaller
  T by its square (y -> y^2, exponents doubled), cost <= 3 units. This depends on the cycle
  only, so merged walks stay merged. A checkpoint T is stored every 64 steps and compared at
  every step; a repeat means a cycle of length <= 64, which is traversed once more to its minimum-T
  class that is then squared (waste <= 200 steps). Every other walk that sees no DP within
  20/theta = 2^64.3 steps is restarted from a fresh uniform pair (stuck-walk rule).
* Total overhead factor on the step count: eps = 5/(12R) + 250/R^2 = 6.5e-4 (charged).

**3.4 Caps and the total-work abort.** The algorithm meters its own charged work with the prices of
section 5 (1.125 units per multiply step, 0.585 per squaring step, the stated costs of escapes,
starts, report handling and verifications) and **aborts with failure as soon as the metered work exceeds
W_cap = 2^125.74**, the submitted bound. The claimed bound therefore holds for every run by construction;
what remains is the probability that the abort fires, which is part of the failure budget below.

*Heuristic H1+ (no drift of the walk statistics).* The cap exceeds the nominal ledger W_nominal = 2^125.7147
by 1.7% (0.0163 units per step; `experiment.py` part D). A multiply step costs 1.125 units and a squaring
step 0.585, so the cap is exceeded only if the fraction of squaring steps among the n <= 2^125.83 steps
falls about 3 percentage points below 3/8 (or the escape events are many times more frequent than the
rates of 3.3). H1+ states that this does not happen: over the whole run the frequencies of the branch bits,
of the 2-cycle events and of the DPs stay within these margins, and the number of reports, false
fingerprint candidates and restarts stays within the caps below, with total probability of any cap or the
abort firing at most 0.01. This is a much weaker statement than independence of the steps; it excludes
only a drift by percentage points over 2^125 steps. Support: (i) the marginal probability of the squaring
branch is exactly 3/8 for every class (T mod 8 is uniform, F2b), so the expected charge of a step is exactly
the nominal one with no independence assumption; (ii) in the toy runs the fluctuation of the mean charge of a
run scales as the independent-steps prediction: rms of (relative deviation x sqrt(steps)) is 0.300, 0.277 and
0.308 observed against 0.283 predicted (300, 150 and 200 runs in q = 67, 163, 211, part B), so there is no excess correlation of
branch labels along the walk at that scale; (iii) the production-field walk of part C has squaring fraction
0.375 +- 0.01. The toy runs cannot test concentration at 2^125 steps; H1+ extrapolates it, exactly as H1
extrapolates the collision law. The success analysis (section 4) charges 0.01 for it.

*Caps.* reports D <= 1.1 S theta + 2^25 (S theta = 2^65.81; 2^25 covers the at most 20 M repeated reports of
the merged trails); at most 2^16 false fingerprint candidates are verified (the expectation is 2^2.7, so
Markov gives probability 2^-13.3; the verification work is at most 2^16 * 600 units = 2^25.2); at most
2^50 restarts (expected about 2^37, Markov 2^-13; each restart is a fresh start of about 2^10 units, 2^60 in
total); the total number of steps <= S(1+eps) + M 20/theta. Repeat steps after the first collision are at
most E = 2^84.32 steps (<= 1.125 units each), negligible against the cap slack 2^119.7. Fruitless cycles
longer than 64 steps need all their steps to be multiply steps (rate below 2^-100 per step, about 2^25.8
events, each wasting at most 20/theta = 2^64.3 steps before the stuck-walk restart: 2^90.1 steps in total).
These caps bound the work and memory stated below.

## 4. Success probability

**H1 (random mapping with branch mix).** The class sequence behaves like iterates of a
random mapping on N classes in which two distinct classes have the same image with
probability (1 - sum_b q_b^2)/N, where q_b are the branch probabilities. Justification: for
c1 != c2 the equation F(c1) = F(c2) has, for each branch pair (b1,b2) with b1 != b2, about one
solution among the N classes (12 candidate representatives, each canonical with probability
1/12); for b1 = b2 = squaring it has none (squaring is a bijection of GT of odd order that
commutes with Frobenius), and for b1 = b2 = a table index it has about 11/12 of one. We
use the conservative value sum q^2 = s^2 + (1-s)^2/R = 0.141006 for all pairs; the exact
value (s^2 + (1-s)^2/(12R)) is smaller by 3.5e-4. Part E of `experiment.py` tests this law directly (section 7). This is the same correction as the Brent-Pollard
factor R/(R-1) of the pure adding walk (s = 0), and for Pollard's original three-branch walk
it is the mechanism that makes Pollard's original square-and-add walk slower than an r-adding walk
(Teske 1998). Hence the effective class count is
N' = N / (1 - sum q^2), log2 N' = 251.4914.

Under H1 the probability that the first S steps contain a class collision is
1 - exp(-S^2/(2N')). Choose P_model = 0.43:
`S = sqrt(2 N' ln(1/0.57)) = 2^125.8302.`
Losses before the first useful collision is detected: detection failure e^-20 = 2e-9, a
collision inside a single walk (fraction about 1/M = 1e-6), beta^1 = beta^2 (2^-255), cap
cap and abort failures (at most 0.01 together, H1+, section 3.4), steps lost to restarts (20 e^-20 = 4.1e-8 of all steps; about 2^37 restarts,
below the cap 2^40, each costing a fresh start of about 2^10 units, 2^47 in total, inside the 2^64.1
bucket of section 5). Their sum is below 0.011. The model value 0.43 therefore leaves 0.03 absolute
slack over the required 0.39 for the inaccuracy of H1. Toy measurements (section 7) of the fraction of
runs with a useful collision within S(0.43) are 0.41-0.475 per setting, pooled 0.447 +- 0.011 over 2,200
runs, at 0.437 for S(0.42); individual settings of a few hundred runs scatter by +-0.04.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Primitive costs, scalar AArch64, using the measured blst_fp_mul = 383 instructions: a
modular add/sub/neg/double of 6 limbs called as a routine, **45** instructions (6 add/sub, 6
conditional-subtract, 6 select, 12 loads, 6 stores, 9 for call, return and address set-up); Fp2
squaring (complex method: (a0+a1)(a0-a1) and 2 a0 a1) = 2 Fp multiplications + 3 additions = 901;
Fp2 addition = 90.

**Table multiply:** the product y*T[h][t] is the dense blst_fp12_mul reference operation,
charged exactly 1.

**Squaring branch:** 9 Fp2 squarings + 34 Fp2 additions (F3) + 200 moves = 9*901 + 34*90 + 200
= 11,369; charged **11,400 instructions = 0.57 units**. A squaring step needs no key computation.

**Per-step work.** Every step (both branches):

| Work | Instructions |
| --- | ---: |
| read Re(c_0) (6 limbs), branch bits, DP test, table address | 50 |
| (alpha, beta) update mod r (two 4-limb add + conditional subtract, or doubling) | 100 |
| 2-cycle and checkpoint comparisons on T, history | 60 |
| loop control | 30 |
| margin for loads and moves | 30 |
| **total counted 270, charged** | **300 = 0.015 units** |

Multiply steps in addition (the keys give kappa):

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a zeta, a Re d_1, b Im d_1, B zeta) | 1,532 |
| 1 Fp addition (B) + 8 Fp subtractions/negations, 45 each | 405 |
| minimum of the 12 keys (12 top-limb loads, 11 x 3 compare/select, tie test, index) | 150 |
| load of the w^1 coefficient, kappa, twist index | 32 |
| **total counted 2,119, charged** | **2,200 = 0.11 units** |

So a multiply step costs 1 + 0.015 + 0.11 = 1.125 units and a squaring step 0.57 + 0.015 = 0.585 units.
Average step cost c = (5/8) * 1.125 + (3/8) * 0.585 = **0.9225 units**.

Setup and output (all <= 2^31 units): pairing g and verification pairing (<= 2^14 each),
check z^r = 1, the table (1024 multipliers x (2 exponentiations <= 1,021 units + 11 Frobenius
maps) <= 2^20.1), M = 2^20 starts g^alpha z^beta (<= 2^30.1), final [x]P0 (<= 510 G1 operations,
each < 1 unit). DP handling: at most 2^65.95 reports (the cap) x 6,000 instructions (the 12 keys for kappa, 2,119; hash
insertion, two 255-bit modular multiplications by p^kappa, memory traffic) / 20,000 = 2^64.2 units. All non-walk
work is therefore below 2^64.2 + 2^31 < 2^64.3 units.

Detection delay: after the first useful collision all M walks continue for 20/theta steps:
E = M * 20 / theta = 2^84.32 steps.

`W = c * S * (1 + eps) + c * E + 2^64.3 = 0.9225 * 2^125.8302 * 1.000645 + 0.9225 * 2^84.32 + 2^64.3 = 2^125.7147`.

Submitted bound: **time_log2 = 125.74** (0.025 bits above the computed value; this covers a 1.7% error
in the average step price). `experiment.py` (part D) recomputes this ledger with exact decimal
arithmetic.

Memory: at most 2^65.95 reports x 80 bytes (128-bit fingerprint, two 256-bit exponents) at load factor
2/3 = 2^72.9 bytes, plus the 7.1 MB table, the 0.8 MB exponent tables and 2^20 walk states (at most
2 KB each, including the remembered previous state): below 2^73 bytes.
Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 above (random-mapping behaviour of the mixed walk on classes), H2 (pseudorandom branch and DP bits
  of T), independence of table indices (used for cycle probabilities). These are standard in
  rho analyses (Teske 1998/2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010), supported here by toy and production-field runs, not proved.
* The squaring and per-step prices follow from the measured fp_mul and a stated primitive model (modular
  add 45 instructions); they are analytic instruction bounds, not measured routines. A pessimistic
  squaring price of 1.0 unit (a full dense multiplication) would remove the gain (the bound would be
  2^126.00, worse than the plain walk); each +0.1 unit on the squaring price adds 0.057 bits, each
  +100 instructions of per-step work 0.007 bits. A Callgrind count of blst's cyclotomic squaring would replace the
  analytic price.
* Collision statistics are validated in toy fields only. The measured mean of the mixed walk is within
  5% of the model in all settings (pooled -1.4%); the plain-walk controls are -3.7% and +4.3%. An independent
  rerun of the script by a separate model instance (q = 67 and 163, about 1,800 mixed-walk runs, several
  seeds, walks and dp_bits) gave a pooled fraction of 0.425 at S(0.42), one 183-run setting 0.344.
* This is a generic square-root attack. It does not use index-calculus (TNFS) methods in Fp12.
  It is a valid upper bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

* **Part A (BLS12-381):** GT membership; Frobenius formula = p-th power; 4-multiplication key
  formula = keys of the 12 explicitly computed conjugates, all distinct; kappa picks one identical
  canonical element from every conjugate; p^6 inverts; ord_r(p) = 12; 1+p^4+p^8 = 0; r-1
  factorisation; ord_r(2); the Granger-Scott formula equals y*y with 9 Fp2 squarings and 34 Fp2
  additions.
* **Part C (BLS12-381, production step rule):** a real walk of thousands of steps in GT with the
  production branch rule (branch and multiplier from the bits of Re(c_0), twisted table, squaring,
  2-cycle escape). Checks: every one of the first 150 steps lands in the class of y^ M_h or
  y^^2 as computed with explicit Frobenius; the exponent invariant y = g^alpha z^beta at the
  end; the squaring fraction (0.375 expected) and chi-square of the low bits of Re(c_0).
* **Part E (pair-collision law, toy fields):** sample m distinct classes, apply one step F to each and
  count coincidences of the images; N times the pair-collision rate is the quantity H1 predicts as
  1 - s^2 - (1-s)^2/(12R). Results (60,000 sampled classes, relative standard error about 0.9% except the last row):

  | q, R, s | observed | exact prediction | conservative value used | random mapping |
  | --- | ---: | ---: | ---: | ---: |
  | 67, 8, 3/8 | 0.8552 | 0.8553 | 0.8105 | 1 |
  | 67, 64, 3/8 | 0.8621 | 0.8589 | 0.8533 | 1 |
  | 67, 8, 0 | 0.9936 | 0.9896 | 0.8750 | 1 |
  | 67, 64, 0 | 1.0013 | 0.9987 | 0.9844 | 1 |
  | 67, 16, 1/2 | 0.7416 | 0.7487 | 0.7344 | 1 |
  | 211, 16, 3/8 (2.1%) | 0.8872 | 0.8573 | 0.8350 | 1 |

  The squaring branch removes exactly its s^2 of the collision probability, and the table branches
  keep 11/12 of theirs (rows with s = 0 and R = 8 separate the exact from the conservative value).
* **Part B (toy fields of the same tower shape, q = 67, 139, 163, 211, 379):** the complete
  algorithm (keys with tie-breaking by full conjugates, kappa trick, twisted table, squaring,
  2-cycle and checkpoint escapes, DPs, collision solving, verification g^x = z). Toy keys take
  few values, so the branch/DP bits are taken from a hash of the canonical conjugate (a
  stand-in for the bits of the 381-bit K; Part C tests those bits in the real field). Part B also reports the
  fluctuation of the mean charge per run against the independent-steps prediction (section 3.4).

Offline aggregate of the same script (all discrete logs recovered and verified; `trials` runs per
setting, `dp_bits` 3 and `walks` 4 unless noted; S(x) is the step budget at model success x):

| q, N classes, R, seed | squaring | runs | mean T_c / sqrt N (model) | P[T_c <= S(0.40)] | P[T_c <= S(0.42)] | P[T_c <= S(0.43)] |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 67, 129,145, 64, seed 11 | 3/8 | 400 | 1.295 (1.357) | 0.438 | 0.468 | 0.475 |
| 211, 721,279, 256, seed 21 | 3/8 | 300 | 1.316 (1.353) | 0.433 | 0.443 | 0.453 |
| 163, 4,524,905, 64, seed 31 | 3/8 | 800 | 1.333 (1.357) | 0.415 | 0.436 | 0.445 |
| 163, 4,524,905, 64, seed 22 | 3/8 | 100 | 1.378 (1.357) | 0.390 | 0.400 | 0.410 |
| 163, 4,524,905, 64, seeds 61-63, walks 8, dp_bits 4 | 3/8 | 600 | 1.380 (1.357) | 0.405 | 0.422 | 0.433 |
| **pooled mixed walk** | 3/8 | 2,200 | 1.338 (1.357) | 0.418 | 0.437 | 0.447 |
| 67, 129,145, 64, seed 12 (control) | 0 | 400 | 1.217 (1.263) | 0.420 | 0.445 | 0.455 |
| 163, 4,524,905, 64, seed 32 (control) | 0 | 600 | 1.318 (1.263) | 0.370 | 0.385 | 0.395 |

The control rows are plain adding walks (their model values use the plain-walk effective class count). The
last row is 1.4 sigma below the 0.413 that the earlier plain-walk record reports for the same field. Measured
2-cycle rate at R = 64 with squaring: 4.5e-4 to 5.4e-4 per step across these settings (predicted
(5/8)^2/(12R) = 5.1e-4, bound 1/(12R) = 1.3e-3); without squaring 1.28e-3 and 1.35e-3 (1.30e-3 predicted).
Reproduce a row with e.g.
`echo '{"q":163,"trials":800,"R":64,"sq8":3,"dp_bits":3,"walks":4,"seed":31,"production_checks":0,"seconds":150}' | python3 -I experiment.py`.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp.
70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Granger & Scott, PKC 2010; Aranha et
al., EUROCRYPT 2011; Bos, Kleinjung & Lenstra, ANTS 2010; mpjunior92, submission 7da5f23a.
