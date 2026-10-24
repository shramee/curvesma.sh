"""Frobenius-class rho with a cyclotomic-squaring branch read from a trace invariant: identities, ledger, toy runs.

Part A  exact identities on BLS12-381 itself (key formula, canonical conjugate,
        Granger-Scott squaring with operation counts, orders of p and 2 modulo r).
Part C  the production step rule (twisted-table multiply, squaring branch, branch and
        distinguished-point bits taken from K) executed in the real BLS12-381 GT for
        a few thousand steps, checked against explicit Frobenius arithmetic.
Part B  the complete algorithm (cycle escapes, distinguished points, collision
        solving) in a small field of the same tower shape, with first-collision
        statistics compared with the random-mapping law adjusted for the branch mix.
Part D  exact recomputation of the work ledger used in proof.md.
Toy runs are evidence at their own scale only; they do not validate production costs.
"""
import json, math, random, sys, time
from fractions import Fraction

T0 = time.time()
cfg = json.load(sys.stdin)
Q_TOY = cfg.get("q", 67)
TRIALS = cfg.get("trials", 20)
R = cfg.get("R", 64)
SQ8 = cfg.get("sq8", 3)
DP_BITS = cfg.get("dp_bits", 3)
WALKS = cfg.get("walks", 4)
SEED = cfg.get("seed", 1)
PROD_CHECKS = cfg.get("production_checks", 2)
PROD_STEPS = cfg.get("production_steps", 3000)
BUDGET = cfg.get("seconds", 12)
PAIR_SAMPLES = cfg.get("pair_samples", 0)
rng = random.Random(SEED)

M64 = (1 << 64) - 1


def mix64(x):
    x &= M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
    return x ^ (x >> 31)


class Tower:
    """F_q2 = F_q[u]/(u^2+1), F_q12 = F_q2[w]/(w^6-xi), xi = 1+u (the blst/profile tower)."""

    def __init__(self, q):
        self.q = q
        self.count = None
        self.d = [self.f2pow((1, 1), m * (q - 1) // 6) for m in range(6)]
        z = self.f2pow((1, 1), (q * q - 1) // 6)
        assert z[1] == 0 and (z[0] * z[0] - z[0] + 1) % q == 0
        self.zeta = z[0]
        self.one = ((1, 0),) + ((0, 0),) * 5
        self._fk = {}

    def f2mul(self, a, b):
        q = self.q
        return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

    def f2pow(self, a, e):
        r = (1, 0)
        while e:
            if e & 1:
                r = self.f2mul(r, a)
            a = self.f2mul(a, a)
            e >>= 1
        return r

    def f2sqr(self, a):
        if self.count is not None:
            self.count["fp2_sqr"] += 1
        q = self.q
        return (((a[0] + a[1]) * (a[0] - a[1])) % q, (2 * a[0] * a[1]) % q)

    def f2add(self, a, b):
        if self.count is not None:
            self.count["fp2_addsub"] += 1
        q = self.q
        return ((a[0] + b[0]) % q, (a[1] + b[1]) % q)

    def f2sub(self, a, b):
        if self.count is not None:
            self.count["fp2_addsub"] += 1
        q = self.q
        return ((a[0] - b[0]) % q, (a[1] - b[1]) % q)

    def f2mulxi(self, a):
        if self.count is not None:
            self.count["fp2_addsub"] += 1
        q = self.q
        return ((a[0] - a[1]) % q, (a[0] + a[1]) % q)

    def mul(self, x, y):
        q = self.q
        acc = [[0, 0] for _ in range(11)]
        for i in range(6):
            a0, a1 = x[i]
            for j in range(6):
                b0, b1 = y[j]
                t = acc[i + j]
                t[0] += a0 * b0 - a1 * b1
                t[1] += a0 * b1 + a1 * b0
        out = []
        for m in range(6):
            c0, c1 = acc[m]
            if m + 6 <= 10:
                h0, h1 = acc[m + 6]
                c0 += h0 - h1
                c1 += h0 + h1
            out.append((c0 % q, c1 % q))
        return tuple(out)

    def pow(self, x, e):
        r = self.one
        while e:
            if e & 1:
                r = self.mul(r, x)
            x = self.mul(x, x)
            e >>= 1
        return r

    def frobk(self, y, k):
        k %= 12
        if k == 0:
            return y
        D = self._fk.get(k)
        if D is None:
            D = [self.f2pow((1, 1), m * (self.q ** k - 1) // 6) for m in range(6)]
            self._fk[k] = D
        q = self.q
        if k & 1:
            return tuple(self.f2mul((c[0], (-c[1]) % q), D[m]) for m, c in enumerate(y))
        return tuple(self.f2mul(c, D[m]) for m, c in enumerate(y))

    def keys(self, y):
        """Re(w^1 coefficient) of y^(q^k), k=0..11, from 4 multiplications in F_q."""
        q, z = self.q, self.zeta
        a = y[1][0]
        b = (y[1][0] * self.d[1][0] + y[1][1] * self.d[1][1]) % q
        za, zb = z * a % q, z * b % q
        A = [a, za, (za - a) % q, (-a) % q, (-za) % q, (a - za) % q]
        B = [b, zb, (zb - b) % q, (-b) % q, (-zb) % q, (b - zb) % q]
        out = [0] * 12
        for j in range(6):
            out[2 * j] = A[j]
            out[2 * j + 1] = B[j]
        return out

    def keys_direct(self, y):
        out = []
        for k in range(12):
            out.append(self.frobk(y, k)[1][0])
        return out

    def fp4sqr(self, a, b):
        a2, b2 = self.f2sqr(a), self.f2sqr(b)
        ab2 = self.f2sub(self.f2sub(self.f2sqr(self.f2add(a, b)), a2), b2)
        return self.f2add(a2, self.f2mulxi(b2)), ab2

    def cyc_sqr(self, y):
        """Granger-Scott squaring, valid on the cyclotomic subgroup (so on GT).
        With A=(c0,c3), B=(c1,c4), C=(c2,c5) in F_q4 = F_q2[v]/(v^2-xi), v=w^3:
          c0'=3A0-2c0, c3'=3A1+2c3, c1'=3(xi C1)+2c1, c4'=3C0-2c4, c2'=3B0-2c2, c5'=3B1+2c5
        where (A0,A1)=A^2, (B0,B1)=B^2, (C0,C1)=C^2."""
        c0, c1, c2, c3, c4, c5 = y
        A0, A1 = self.fp4sqr(c0, c3)
        B0, B1 = self.fp4sqr(c1, c4)
        C0, C1 = self.fp4sqr(c2, c5)
        xC1 = self.f2mulxi(C1)

        def minus(t, u):  # 3t - 2u = t + 2(t-u)
            d = self.f2sub(t, u)
            return self.f2add(t, self.f2add(d, d))

        def plus(t, u):  # 3t + 2u = t + 2(t+u)
            s = self.f2add(t, u)
            return self.f2add(t, self.f2add(s, s))

        return (minus(A0, c0), plus(xC1, c1), minus(B0, c2), plus(A1, c3), minus(C0, c4), plus(B1, c5))

    def rand_gt(self, rng, rp):
        hc = (self.q ** 12 - 1) // rp
        while True:
            x = tuple((rng.randrange(self.q), rng.randrange(self.q)) for _ in range(6))
            y = self.pow(x, hc)
            if y != self.one:
                return y


# ---------------------------------------------------------------- production constants
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
R_FACT = {2: 32, 3: 1, 11: 1, 19: 1, 10177: 1, 125527: 1, 859267: 1, 906349: 2,
          2508409: 1, 2529403: 1, 52437899: 1, 254760293: 2}


def mult_order(a, n, fact):
    o = n
    for pf, e in fact.items():
        for _ in range(e):
            if pow(a, o // pf, n) == 1:
                o //= pf
            else:
                break
    return o


def part_A(rng, checks):
    P = Tower(p)
    lam = p % r
    prod_fact = 1
    for k, e in R_FACT.items():
        prod_fact *= k ** e
    o2 = mult_order(2, r, {k: e for k, e in R_FACT.items()})
    ol = mult_order(lam, r, R_FACT)
    t_min = o2 // math.gcd(o2, 12)
    out = {
        "r_minus_1_factorisation_verified": prod_fact == r - 1,
        "order_of_p_mod_r": ol,
        "lambda6_is_minus1": pow(lam, 6, r) == r - 1,
        "one_plus_l4_plus_l8_is_0": (1 + pow(lam, 4, r) + pow(lam, 8, r)) % r == 0,
        "log2_order_of_2_mod_r": round(math.log2(o2), 3),
        "smallest_t_with_2^t_in_<p>_log2": round(math.log2(t_min), 3),
        "zeta_primitive_6th_root": (P.zeta * P.zeta - P.zeta + 1) % p == 0,
    }
    ok = dict.fromkeys(["gt", "frob", "keys", "distinct", "canon", "inv", "sq", "trace_inv", "trace_diff"], True)
    traces = []
    counts = []
    hcof = (p ** 12 - 1) // r
    for _ in range(checks):
        x = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
        y = P.pow(x, hcof)
        ok["gt"] &= P.pow(y, r) == P.one and y != P.one
        ok["frob"] &= P.frobk(y, 1) == P.pow(y, p)
        kd, kf = P.keys_direct(y), P.keys(y)
        ok["keys"] &= kd == kf
        ok["distinct"] &= len(set(kf)) == 12
        conj = [P.frobk(y, k) for k in range(12)]
        canon = set()
        for c in conj:
            kk = P.keys(c)
            kap = min(range(12), key=lambda i: kk[i])
            canon.add(P.frobk(c, kap))
        ok["canon"] &= len(canon) == 1
        ok["inv"] &= P.mul(conj[6], y) == P.one
        ok["trace_inv"] &= len({c[0][0] for c in conj}) == 1       # Re(c_0) is the same for all 12 conjugates
        traces.append(y[0][0])
        P.count = {"fp2_sqr": 0, "fp2_addsub": 0}
        s = P.cyc_sqr(y)
        counts.append(dict(P.count))
        P.count = None
        ok["sq"] &= s == P.mul(y, y)
    ok["trace_diff"] &= len(set(traces)) == len(traces)
    out.update({"random_GT_elements": checks, **{k + "_ok": v for k, v in ok.items()}})
    out["granger_scott_op_counts"] = counts[0] if counts else None
    return out


def part_C(rng, steps, rprod=16, verify=150):
    """Production step rule in the real GT: twisted-table multiply / squaring, bits from T = Re(c_0)."""
    P = Tower(p)
    hcof = (p ** 12 - 1) // r
    lam = [pow(p, k, r) for k in range(12)]
    g = P.pow(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6)), hcof)
    x0 = rng.randrange(1, 1 << 40)
    z0 = P.pow(g, x0)
    ah = [rng.randrange(1, 1 << 40) for _ in range(rprod)]
    bh = [rng.randrange(1, 1 << 40) for _ in range(rprod)]
    Mh = [P.mul(P.pow(g, ah[i]), P.pow(z0, bh[i])) for i in range(rprod)]
    Tab = [[P.frobk(Mh[i], k) for k in range(12)] for i in range(rprod)]
    # exponents are tracked exactly: the start element is g^a0 z0^b0
    a0, b0 = rng.randrange(1, r), rng.randrange(1, r)
    y = P.mul(P.pow(g, a0), P.pow(z0, b0))
    al, be = a0, b0
    n_sq = n_mul = bad_class = 0
    bins_sq = [0] * 8
    bins_h = [0] * 16
    two_cycles = 0
    hist = []
    for t in range(steps):
        u = y[0][0]                                      # class invariant T = Re(c_0): no multiplication needed
        if len(hist) >= 2 and hist[-2][0] == u:          # fruitless 2-cycle: square the smaller-T class
            two_cycles += 1
            _, yy, aa, bb = min(hist[-1], (u, y, al, be), key=lambda s: s[0])
            y, al, be = P.cyc_sqr(yy), 2 * aa % r, 2 * bb % r
            hist = []
            continue
        hist.append((u, y, al, be))
        if len(hist) > 2:
            hist.pop(0)
        br = u & 7
        bins_sq[br] += 1
        if br < SQ8:
            y2 = P.cyc_sqr(y)                            # squaring step: no keys, no twist
            if t < verify:
                kap = min(range(12), key=lambda i: P.keys(y)[i])
                bad_class += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), P.frobk(y, kap))
            y, al, be = y2, 2 * al % r, 2 * be % r
            n_sq += 1
        else:
            kf = P.keys(y)                               # multiply step: 4 Fp multiplications give kappa
            kap = kf.index(min(kf))
            h = (u >> 3) % rprod
            bins_h[(u >> 3) & 15] += 1
            idx = (12 - kap) % 12
            y2 = P.mul(y, Tab[h][idx])
            if t < verify:
                bad_class += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), Mh[h])
            y = y2
            al = (al + ah[h] * lam[idx]) % r
            be = (be + bh[h] * lam[idx]) % r
            n_mul += 1
    exp_ok = y == P.mul(P.pow(g, al), P.pow(z0, be))
    exp_sq = sum(bins_sq[:SQ8])
    n_br = sum(bins_h)
    n_all = sum(bins_sq)
    chi_sq = sum((c - n_all / 8) ** 2 / (n_all / 8) for c in bins_sq)
    chi_h = sum((c - n_br / 16) ** 2 / (n_br / 16) for c in bins_h) if n_br else None
    return {"steps": steps, "R_table": rprod, "squaring_steps": n_sq, "multiply_steps": n_mul,
            "squaring_fraction": round(n_sq / max(1, n_all), 4), "target_fraction": SQ8 / 8,
            "explicit_class_check_failures_in_first_%d" % verify: int(bad_class),
            "exponent_invariant_y_equals_g^al_z^be": exp_ok,
            "chi2_low3_bits_of_Re_c0_7dof": round(chi_sq, 2),
            "chi2_next4_bits_of_Re_c0_15dof": round(chi_h, 2) if chi_h is not None else None,
            "two_cycle_events": two_cycles}


# ---------------------------------------------------------------- toy algorithm (Part B)
def toy_setup(q):
    phi = q ** 4 - q * q + 1
    rp = phi
    d = 2
    while d * d <= rp:
        while rp % d == 0 and rp > d:
            rp //= d
        d += 1
    return rp


def classinfo(T, y):
    """kappa, canonical conjugate, minimum key. Production never builds the canonical
    element; the toy does, and breaks equal-key ties by comparing full conjugates."""
    kf = T.keys(y)
    m = min(kf)
    ks = [k for k in range(12) if kf[k] == m]
    if len(ks) == 1:
        kap = ks[0]
        can = T.frobk(y, kap)
    else:
        can, kap = min((T.frobk(y, k), k) for k in ks)
    return kap, can, m


def toy_hash(can):
    h = 0x9E3779B97F4A7C15
    for c in can:
        h = mix64(h ^ (c[0] * 0x100000001B3 + c[1]))
    return h


def run_toy(T, q, rp, rng, walks, Rm, sq8, dp_bits, cap_factor=40):
    qk = [pow(q, k, rp) for k in range(12)]
    hc = (q ** 12 - 1) // rp
    g = T.rand_gt(rng, rp)
    xs = rng.randrange(1, rp)
    z0 = T.pow(g, xs)
    ah = [rng.randrange(rp) for _ in range(Rm)]
    bh = [rng.randrange(rp) for _ in range(Rm)]
    Mh = [T.mul(T.pow(g, ah[i]), T.pow(z0, bh[i])) for i in range(Rm)]
    Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
    At = [[ah[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    Bt = [[bh[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    mask = (1 << dp_bits) - 1

    def new_walk():
        a, b = rng.randrange(rp), rng.randrange(rp)
        return {"y": T.mul(T.pow(g, a), T.pow(z0, b)), "a": a, "b": b, "hist": [], "ck": None, "n": 0, "since": 0}

    def advance(y, a, b):
        kap, can, K = classinfo(T, y)
        u = toy_hash(can)
        if (u & 7) < sq8:
            return T.cyc_sqr(y), 2 * a % rp, 2 * b % rp, kap, can, K, u
        h = (u >> 3) % Rm
        idx = (12 - kap) % 12
        return T.mul(y, Tab[h][idx]), (a + At[h][idx]) % rp, (b + Bt[h][idx]) % rp, kap, can, K, u

    W = [new_walk() for _ in range(walks)]
    visited, dps = {}, {}
    steps, Tc, found = 0, None, None
    stats = {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0, "sq": 0, "ord": 0}
    while found is None:
        for w in W:
            y, a, b = w["y"], w["a"], w["b"]
            kap, can, K = classinfo(T, y)
            u = toy_hash(can)
            ca, cb = a * qk[kap] % rp, b * qk[kap] % rp
            if Tc is None:
                old = visited.get(can)
                if old is None:
                    visited[can] = (ca, cb)
                elif old != (ca, cb):
                    Tc = steps
            if (u >> 13) & mask == 0:
                old = dps.get(can)
                if old is not None and old[1] != cb:
                    xr = (old[0] - ca) * pow(cb - old[1], -1, rp) % rp
                    if T.pow(g, xr) != z0:
                        raise RuntimeError('recovered logarithm failed verification')
                    found = xr
                    break
                if old is not None:
                    stats["useless"] += 1
                dps[can] = (ca, cb)
                w["since"] = 0
            hist = w["hist"]
            key = (K, can)
            if len(hist) >= 2 and hist[-2][0][1] == can:
                stats["two_cycles"] += 1
                s1, s2 = hist[-1], (key, y, a, b)
                _, yy, aa, bb = min(s1, s2, key=lambda s: s[0])
                w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
                w["hist"], w["ck"] = [], None
                steps += 1
                continue
            if w["ck"] is not None and w["ck"][0] == can and w["ck"][1] < w["n"]:
                stats["long_cycles"] += 1
                best = (key, y, a, b)
                yy, aa, bb = y, a, b
                for _ in range(w["n"] - w["ck"][1]):
                    yy, aa, bb, kp, cn, Kn, un = advance(yy, aa, bb)
                    steps += 1
                    k2, c2, K2 = classinfo(T, yy)
                    if (K2, c2) < best[0]:
                        best = ((K2, c2), yy, aa, bb)
                _, yy, aa, bb = best
                w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
                w["hist"], w["ck"] = [], None
                steps += 1
                continue
            if w["n"] % 64 == 0:
                w["ck"] = (can, w["n"])
            hist.append((key, y, a, b))
            if len(hist) > 2:
                hist.pop(0)
            stats["ord"] += 1
            stats["sq"] += (u & 7) < sq8
            w["y"], w["a"], w["b"], _, _, _, _ = advance(y, a, b)
            w["n"] += 1
            w["since"] += 1
            steps += 1
            if w["since"] > cap_factor << dp_bits:
                stats["restarts"] += 1
                w.update(new_walk())
    return Tc, steps, stats


def part_B(rng, q, trials, Rm, sq8, dp_bits, walks, budget):
    T = Tower(q)
    rp = toy_setup(q)
    if not (pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))):
        raise RuntimeError('toy field does not have embedding degree 12')
    N = (rp - 1) // 12
    s = sq8 / 8
    sumq2 = s * s + (1 - s) ** 2 / Rm
    Np = N / (1 - sumq2)
    Tcs, tot, ev, dev = [], [], {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0}, []
    t_start = time.time()
    while len(Tcs) < trials and time.time() - t_start < budget:
        Tc, steps, st = run_toy(T, q, rp, rng, walks, Rm, sq8, dp_bits)
        Tcs.append(Tc)
        tot.append(steps)
        for k in ev:
            ev[k] += st[k]
        sf = sq8 / 8
        nominal = (1 - sf) * 1.125 + sf * 0.585
        charge = ((st["ord"] - st["sq"]) * 1.125 + st["sq"] * 0.585) / max(1, st["ord"])
        dev.append((charge / nominal - 1, st["ord"], nominal))
    n = len(Tcs)
    sN = math.sqrt(N)
    pred = math.sqrt(math.pi * Np / 2)
    s40 = math.sqrt(2 * Np * math.log(1 / 0.60))
    s42 = math.sqrt(2 * Np * math.log(1 / 0.58))
    s43 = math.sqrt(2 * Np * math.log(1 / 0.57))
    allsteps = max(1, sum(tot))
    return {"q": q, "r_prime": rp, "classes_N": N, "R": Rm, "sq8": sq8, "squaring_prob": s,
            "sum_q2": round(sumq2, 6), "effective_N_over_N": round(1 / (1 - sumq2), 4),
            "trials_completed": n, "logs_recovered_and_verified": n,
            "mean_Tc_over_sqrtN": round(sum(Tcs) / n / sN, 4) if n else None,
            "model_mean_over_sqrtN": round(pred / sN, 4),
            "fraction_Tc_le_S_model_0.40": round(sum(t <= s40 for t in Tcs) / max(1, n), 3),
            "fraction_Tc_le_S_model_0.42": round(sum(t <= s42 for t in Tcs) / max(1, n), 3),
            "fraction_Tc_le_S_model_0.43": round(sum(t <= s43 for t in Tcs) / max(1, n), 3),
            "S_model_0.40_over_sqrtN": round(s40 / sN, 4), "S_model_0.42_over_sqrtN": round(s42 / sN, 4),
            "S_model_0.43_over_sqrtN": round(s43 / sN, 4),
            "mean_total_steps_incl_DP_delay_over_sqrtN": round(sum(tot) / max(1, n) / sN, 4),
            "charge_fluctuation_rms_of_(dev*sqrt(steps))_observed": round(math.sqrt(sum(d * d * k for d, k, _ in dev) / max(1, len(dev))), 4),
            "charge_fluctuation_iid_prediction": round(math.sqrt(s * (1 - s)) * 0.54 / dev[0][2], 4) if dev else None,
            "events": ev, "two_cycles_per_step": ev["two_cycles"] / allsteps,
            "Tc_over_sqrtN": [round(t / sN, 4) for t in Tcs] if cfg.get("dump") else None}


# ---------------------------------------------------------------- pair-collision law (Part E)
def part_E(rng, q, m, Rm, sq8):
    """Sample m distinct classes, apply one step F to each, and count image coincidences.
    H1 predicts N * P[F(c1) = F(c2)] = 1 - s^2 - (1-s)^2/(12R) for distinct classes c1, c2."""
    T = Tower(q)
    rp = toy_setup(q)
    N = (rp - 1) // 12
    g = T.rand_gt(rng, rp)
    Mh = [T.pow(g, rng.randrange(1, rp)) for _ in range(Rm)]
    Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
    gp = [g]
    for _ in range(rp.bit_length()):
        gp.append(T.mul(gp[-1], gp[-1]))
    classes = {}
    while len(classes) < m:
        e = rng.randrange(1, rp)
        y, j = T.one, 0
        while e:
            if e & 1:
                y = T.mul(y, gp[j])
            e >>= 1
            j += 1
        kap, can, _ = classinfo(T, y)
        if can not in classes:
            classes[can] = (y, kap)
    images = {}
    for can, (y, kap) in classes.items():
        u = toy_hash(can)
        if (u & 7) < sq8:
            y2 = T.cyc_sqr(y)
        else:
            y2 = T.mul(y, Tab[(u >> 3) % Rm][(12 - kap) % 12])
        _, can2, _ = classinfo(T, y2)
        images[can2] = images.get(can2, 0) + 1
    pairs = sum(n * (n - 1) // 2 for n in images.values())
    tot = m * (m - 1) // 2
    s = sq8 / 8
    obs = pairs / tot * N
    return {"q": q, "classes_N": N, "sampled_classes": m, "R": Rm, "squaring_prob": s,
            "image_coincidence_pairs": pairs,
            "N_times_pair_collision_rate_observed": round(obs, 4),
            "relative_standard_error": round(1 / math.sqrt(max(1, pairs)), 4),
            "predicted_exact_1_minus_s2_minus_(1-s)^2/(12R)": round(1 - s * s - (1 - s) ** 2 / (12 * Rm), 4),
            "predicted_conservative_1_minus_s2_minus_(1-s)^2/R": round(1 - s * s - (1 - s) ** 2 / Rm, 4),
            "random_mapping": 1.0}


# ---------------------------------------------------------------- ledger (Part D)
def ledger():
    """Exact recomputation of the work ledger of proof.md section 5 (instruction counts / C, C = 20,000)."""
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    ln2 = Decimal(2).ln()
    lg = lambda v: float(Decimal(v).ln() / ln2)
    C = 20000
    FPMUL, FPADD = 383, 45                      # blst_fp_mul (measured); modular add/sub/neg as a called routine
    FP2SQR = 2 * FPMUL + 3 * FPADD              # complex squaring: 2 Fp mults + 3 Fp add/sub
    FP2ADD = 2 * FPADD
    GS = 9 * FP2SQR + 34 * FP2ADD + 200         # Granger-Scott: 9 Fp2 sqr + 34 Fp2 add/sub (Part A) + moves
    common_items = {"read Re(c_0), branch bits, DP test, table address": 50,
                    "(alpha, beta) update mod r (add or doubling)": 100,
                    "2-cycle and checkpoint comparisons on T, history": 60,
                    "loop control": 30,
                    "margin for loads and moves": 30}
    mult_items = {"4 Fp multiplications (a*zeta, a*d0, b*d1, B*zeta)": 4 * FPMUL,
                  "1 Fp addition (B) + 8 Fp sub/neg": 9 * FPADD,
                  "minimum of 12 keys (top limbs, tie test)": 150,
                  "load c_1, kappa, twist index": 32}
    common, extra = sum(common_items.values()), sum(mult_items.values())
    SQ_CHARGE, COMMON_CHARGE, MULT_EXTRA_CHARGE = 11400, 300, 2200
    if not (GS <= SQ_CHARGE and common <= COMMON_CHARGE and extra <= MULT_EXTRA_CHARGE):
        raise RuntimeError('charged prices below derived counts')
    N = Fraction(r - 1, 12)
    s, Rt = Fraction(3, 8), 1024
    sumq2 = s * s + (1 - s) ** 2 / Rt
    Np = N / (1 - sumq2)
    Pm = Decimal("0.43")
    Npd = Decimal(Np.numerator) / Decimal(Np.denominator)
    S = (2 * Npd * (1 / (1 - Pm)).ln()).sqrt()
    cbar = (1 - s) * (1 + Fraction(COMMON_CHARGE + MULT_EXTRA_CHARGE, C)) + s * (Fraction(SQ_CHARGE + COMMON_CHARGE, C))
    cb = Decimal(cbar.numerator) / Decimal(cbar.denominator)
    eps = Decimal(5) / Decimal(12 * Rt) + Decimal(250) / Decimal(Rt * Rt)
    E = Decimal(2 ** 20) * 20 * Decimal(2 ** 60)
    SETUP = Decimal(2) ** Decimal("64.3")      # non-walk work, proof.md section 5
    W = cb * S * (1 + eps) + cb * E + SETUP
    # total-work abort (proof.md 3.4): plain arithmetic on the slack between the cap and the nominal ledger
    Wcap = Decimal(2) ** Decimal("125.74")
    n = S * (1 + eps)
    slack = (Wcap - cb * E - SETUP) / n - cb                 # per-step slack of the cap over the nominal ledger
    tol_sq_fraction = slack / Decimal("0.54")                 # tolerated drop of the squaring fraction below 3/8 (charges 1.125 vs 0.585)
    return {"fp2_sqr_instr": FP2SQR, "cyclotomic_square_instr_derived": GS, "cyclotomic_square_charged": SQ_CHARGE,
            "common_items": common_items, "multiply_only_items": mult_items, "common_counted": common, "common_charged": COMMON_CHARGE,
            "multiply_extra_counted": extra, "multiply_extra_charged": MULT_EXTRA_CHARGE,
            "sum_q2": float(sumq2), "log2_N": lg(Decimal(N.numerator) / Decimal(N.denominator)),
            "log2_N_effective": lg(Npd), "log2_S_model_0.43": lg(S), "cbar_units_per_step": float(cb),
            "fruitless_cycle_eps": float(eps), "log2_detection_window_steps": lg(E), "log2_W": lg(W),
            "work_abort_slack_per_step": float(slack), "tolerated_squaring_fraction_drop": float(tol_sq_fraction)}


out = {}
out["partA_production_BLS12_381"] = part_A(rng, PROD_CHECKS)
if cfg.get("production_steps", 3000) > 0 and PROD_CHECKS > 0:
    out["partC_production_walk"] = part_C(rng, PROD_STEPS)
out["partB_toy"] = part_B(rng, Q_TOY, TRIALS, R, SQ8, DP_BITS, WALKS, BUDGET)
if PAIR_SAMPLES > 0:
    out["partE_pair_collision_law"] = part_E(rng, Q_TOY, PAIR_SAMPLES, R, SQ8)
out["partD_ledger"] = ledger()
out["seconds"] = round(time.time() - T0, 2)
out["scope"] = ("Parts A and C check exact identities and the production step rule on real BLS12-381 "
                "GT elements; Part B is a toy field; none validates production costs.")
print(json.dumps(out, indent=1))
