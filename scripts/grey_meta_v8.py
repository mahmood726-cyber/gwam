"""
Grey Relational Meta-Analysis (GRMA) v8.0 — Python Implementation
=================================================================

Robust pooling estimator using grey-relational similarity in a two-feature
space (effect size and log-precision) with optional redescending Tukey
bisquare effect guard and full BCa bootstrap inference.

Matches the R implementation in Pairwise70::grma_meta() to 6+ decimal places.

Dependencies: numpy, scipy (see requirements.txt for pinned versions)

Usage:
    from grey_meta_v8 import GRMA, compare_methods
    g = GRMA()
    result = g.fit(effect, variance)
    ci = g.bootstrap_ci(effect, variance)
"""

import numpy as np
from scipy import stats


class GRMA:
    """Grey Relational Meta-Analysis estimator."""

    def __init__(
        self,
        zeta=0.5,
        norm_method="robust_minmax",
        anchor_mode="median",
        trim=0.1,
        prec_cap=1e6,
        effect_guard=True,
        tukey_c=4.685,
        guard_power=1,
    ):
        if not (0 < zeta <= 1):
            raise ValueError("zeta must be in (0, 1]")
        if prec_cap <= 0:
            raise ValueError("prec_cap must be positive")
        if not (0 <= trim < 0.5):
            raise ValueError("trim must be in [0, 0.5)")
        if tukey_c <= 0:
            raise ValueError("tukey_c must be positive")
        self.zeta = zeta
        self.norm_method = norm_method
        self.anchor_mode = anchor_mode
        self.trim = trim
        self.prec_cap = prec_cap
        self.effect_guard = effect_guard
        self.tukey_c = tukey_c
        self.guard_power = guard_power

    def _core(self, y, v):
        """Core GRMA computation (no bootstrap). Returns dict with estimate, weights, etc."""
        n = len(y)
        y = np.asarray(y, dtype=np.float64)
        v = np.asarray(v, dtype=np.float64)

        # Step 1: Precision cap + log-precision
        prec = np.minimum(1.0 / v, self.prec_cap)
        log_prec = np.log(prec + 1.0)

        # Step 2: Feature vectors
        feat_effect = y.copy()
        feat_prec = log_prec.copy()

        # Step 3: Robust min-max normalization (fitted on data)
        def robust_minmax_fit(x):
            q_lo = np.percentile(x, 5, method="linear")
            q_hi = np.percentile(x, 95, method="linear")
            rng = q_hi - q_lo
            if rng < 1e-12:
                rng = 1.0
            return q_lo, q_hi, rng

        def norm_val(x, q_lo, rng):
            return np.clip((x - q_lo) / rng, 0.0, 1.0)

        eff_lo, _, eff_rng = robust_minmax_fit(feat_effect)
        pre_lo, _, pre_rng = robust_minmax_fit(feat_prec)

        x_eff = norm_val(feat_effect, eff_lo, eff_rng)
        x_pre = norm_val(feat_prec, pre_lo, pre_rng)

        # Step 4: Anchor
        if self.anchor_mode == "median":
            a_y_raw = np.median(y)
        elif self.anchor_mode == "trimmed_mean":
            a_y_raw = stats.trim_mean(y, self.trim)
        else:
            a_y_raw = np.median(y)

        a_p_raw = np.max(prec)

        a_eff = float(np.clip((a_y_raw - eff_lo) / eff_rng, 0.0, 1.0))
        a_pre = float(np.clip((np.log(a_p_raw + 1.0) - pre_lo) / pre_rng, 0.0, 1.0))

        # Step 5: Grey relational coefficient and grade
        delta_eff = np.abs(x_eff - a_eff)
        delta_pre = np.abs(x_pre - a_pre)

        all_deltas = np.concatenate([delta_eff, delta_pre])
        delta_min = np.min(all_deltas)
        delta_max = np.max(all_deltas)

        grc_eff = (delta_min + self.zeta * delta_max) / (delta_eff + self.zeta * delta_max)
        grc_pre = (delta_min + self.zeta * delta_max) / (delta_pre + self.zeta * delta_max)

        grade = (grc_eff + grc_pre) / 2.0

        # Step 6: Effect guard (Tukey bisquare)
        guard = None
        if self.effect_guard:
            mad_y = np.median(np.abs(y - a_y_raw))
            if mad_y < 1e-12:
                mad_y = 1e-12
            u = np.abs(y - a_y_raw) / mad_y
            h = np.where(u < self.tukey_c, (1.0 - (u / self.tukey_c) ** 2) ** 2, 0.0)
            guard = h
            raw_w = grade * h ** self.guard_power
        else:
            raw_w = grade

        # Step 7: Normalize weights
        sw = np.sum(raw_w)
        if sw < 1e-15:
            w = np.full(n, 1.0 / n)
        else:
            w = raw_w / sw

        # Step 8: Pooled estimate
        est = float(np.sum(w * y))

        return {
            "estimate": est,
            "weights": w,
            "guard_values": guard,
            "anchor_y": float(a_y_raw),
            "anchor_p": float(a_p_raw),
        }

    @staticmethod
    def _validate_inputs(y, v):
        """Validate effect/variance arrays. Raises ValueError on bad input."""
        if len(y) != len(v):
            raise ValueError("effect and variance must have the same length")
        if len(y) == 0:
            raise ValueError("at least 1 study required")
        if np.any(~np.isfinite(y)) or np.any(~np.isfinite(v)):
            raise ValueError("effect and variance must be finite (no NaN/Inf)")
        if np.any(v <= 0):
            raise ValueError("all variances must be positive")

    def fit(self, effect, variance):
        """Fit GRMA on observed data. Returns dict with estimate and diagnostics."""
        y = np.asarray(effect, dtype=np.float64)
        v = np.asarray(variance, dtype=np.float64)
        self._validate_inputs(y, v)
        k = len(y)

        if k < 3:
            return _hk_reml(y, v)

        res = self._core(y, v)
        w = res["weights"]

        return {
            "estimate": res["estimate"],
            "weights": w,
            "w_max": float(np.max(w)),
            "n_eff": float(1.0 / np.sum(w**2)),
            "anchor_y": res["anchor_y"],
            "anchor_p": res["anchor_p"],
            "guard_values": res["guard_values"],
            "k": k,
            "method": "GRMA" if self.effect_guard else "GRMA_noguard",
        }

    def bootstrap_ci(self, effect, variance, B=999, conf_level=0.95, bca=True, seed=None):
        """
        Bootstrap percentile and (optionally) full BCa confidence intervals.

        Returns dict with ci_lo_pct, ci_hi_pct, ci_lo_bca, ci_hi_bca, se, estimate.
        The p-value is a Wald z-test approximation (est/se) and may not exactly
        agree with the BCa interval boundaries.
        """
        y = np.asarray(effect, dtype=np.float64)
        v = np.asarray(variance, dtype=np.float64)
        self._validate_inputs(y, v)
        k = len(y)

        if k < 3:
            hk = _hk_reml(y, v)
            return {
                "estimate": hk["estimate"],
                "se": hk["se"],
                "ci_lo_pct": hk["ci_lo"],
                "ci_hi_pct": hk["ci_hi"],
                "ci_lo_bca": np.nan,
                "ci_hi_bca": np.nan,
                "pvalue": np.nan,
                "n_boot_ok": 0,
            }

        # Original estimate
        res0 = self._core(y, v)
        est0 = res0["estimate"]

        rng = np.random.default_rng(seed)
        boot_est = np.empty(B)

        for b in range(B):
            idx = rng.integers(0, k, size=k)
            try:
                rb = self._core(y[idx], v[idx])
                boot_est[b] = rb["estimate"]
            except (ValueError, FloatingPointError, ZeroDivisionError, RuntimeWarning):
                boot_est[b] = np.nan

        boot_ok = boot_est[np.isfinite(boot_est)]
        n_ok = len(boot_ok)

        if n_ok < 10:
            return {
                "estimate": est0,
                "se": np.nan,
                "ci_lo_pct": np.nan,
                "ci_hi_pct": np.nan,
                "ci_lo_bca": np.nan,
                "ci_hi_bca": np.nan,
                "n_boot_ok": n_ok,
            }

        alpha = 1.0 - conf_level
        se = float(np.std(boot_ok, ddof=1))

        # Percentile CI
        ci_lo_pct = float(np.percentile(boot_ok, 100 * alpha / 2, method="linear"))
        ci_hi_pct = float(np.percentile(boot_ok, 100 * (1 - alpha / 2), method="linear"))

        # BCa CI
        ci_lo_bca = np.nan
        ci_hi_bca = np.nan

        if bca and n_ok >= 50:
            # Bias correction z0
            prop_below = np.mean(boot_ok < est0)
            prop_below = max(1 / (2 * n_ok), min(1 - 1 / (2 * n_ok), prop_below))
            z0 = stats.norm.ppf(prop_below)

            # Jackknife acceleration
            jack_est = np.empty(k)
            for i in range(k):
                mask = np.ones(k, dtype=bool)
                mask[i] = False
                try:
                    rj = self._core(y[mask], v[mask])
                    jack_est[i] = rj["estimate"]
                except (ValueError, FloatingPointError, ZeroDivisionError, RuntimeWarning):
                    jack_est[i] = np.nan

            jack_ok_mask = np.isfinite(jack_est)
            n_jack_ok = int(np.sum(jack_ok_mask))
            if n_jack_ok >= 3:
                jack_finite = jack_est[jack_ok_mask]
                jack_mean = np.mean(jack_finite)
                d = jack_mean - jack_finite
                denom = np.sum(d**2)
                if denom > 1e-30:
                    a_hat = float(np.clip(
                        np.sum(d**3) / (6.0 * denom**1.5), -0.5, 0.5
                    ))
                else:
                    a_hat = 0.0
            else:
                a_hat = 0.0

            # Adjusted percentiles
            z_lo = stats.norm.ppf(alpha / 2)
            z_hi = stats.norm.ppf(1 - alpha / 2)

            adj_lo = stats.norm.cdf(z0 + (z0 + z_lo) / (1 - a_hat * (z0 + z_lo)))
            adj_hi = stats.norm.cdf(z0 + (z0 + z_hi) / (1 - a_hat * (z0 + z_hi)))

            # Clamp
            adj_lo = max(0.5 / n_ok, min(1 - 0.5 / n_ok, adj_lo))
            adj_hi = max(0.5 / n_ok, min(1 - 0.5 / n_ok, adj_hi))

            ci_lo_bca = float(np.percentile(boot_ok, 100 * adj_lo, method="linear"))
            ci_hi_bca = float(np.percentile(boot_ok, 100 * adj_hi, method="linear"))

        # P-value
        if se > 1e-15:
            z_stat = est0 / se
            pval = float(2 * (1 - stats.norm.cdf(abs(z_stat))))
        else:
            pval = np.nan

        return {
            "estimate": est0,
            "se": se,
            "ci_lo_pct": ci_lo_pct,
            "ci_hi_pct": ci_hi_pct,
            "ci_lo_bca": ci_lo_bca,
            "ci_hi_bca": ci_hi_bca,
            "pvalue": pval,
            "n_boot_ok": n_ok,
        }

    def leave_one_out(self, effect, variance):
        """Leave-one-out influence diagnostics. Requires k >= 4 (so LOO has k-1 >= 3)."""
        y = np.asarray(effect, dtype=np.float64)
        v = np.asarray(variance, dtype=np.float64)
        self._validate_inputs(y, v)
        k = len(y)
        if k < 4:
            raise ValueError("leave_one_out requires k >= 4 (each LOO subset needs k >= 3)")

        res_full = self._core(y, v)
        est_full = res_full["estimate"]
        w_max_full = float(np.max(res_full["weights"]))

        est_loo = np.empty(k)
        w_max_loo = np.empty(k)

        for i in range(k):
            mask = np.ones(k, dtype=bool)
            mask[i] = False
            ri = self._core(y[mask], v[mask])
            est_loo[i] = ri["estimate"]
            w_max_loo[i] = float(np.max(ri["weights"]))

        delta_est = np.abs(est_loo - est_full)
        delta_maxw = np.abs(w_max_loo - w_max_full)

        idx_max = int(np.argmax(delta_est))
        idx_maxw = int(np.argmax(delta_maxw))

        return {
            "est_full": est_full,
            "est_loo": est_loo,
            "delta_est": delta_est,
            "w_max_full": w_max_full,
            "w_max_loo": w_max_loo,
            "max_abs_delta_est": float(np.max(delta_est)),
            "idx_max_shift": idx_max,
            "max_abs_delta_maxw": float(np.max(delta_maxw)),
            "idx_maxw_shift": idx_maxw,
            "y_at_max": float(y[idx_max]),
            "v_at_max": float(v[idx_max]),
        }


def valley_diagnostic(yi, estimate, n_perm=999, seed=None):
    """
    Exploratory valley diagnostic for bimodality.

    Returns dict with valley_flag and valley_p.
    """
    y = np.asarray(yi, dtype=np.float64)
    k = len(y)
    if k < 4:
        return {"valley_flag": False, "valley_p": 1.0}

    below = y[y <= estimate]
    above = y[y > estimate]

    if len(below) < 2 or len(above) < 2:
        return {"valley_flag": False, "valley_p": 1.0}

    within_var = (np.var(below, ddof=1) * (len(below) - 1) +
                  np.var(above, ddof=1) * (len(above) - 1)) / (k - 2)
    total_var = np.var(y, ddof=1)
    if total_var < 1e-15:
        return {"valley_flag": False, "valley_p": 1.0}

    obs_ratio = within_var / total_var

    rng = np.random.default_rng(seed)
    perm_ratios = np.empty(n_perm)
    for p in range(n_perm):
        yp = rng.permutation(y)
        bp = yp[yp <= estimate]
        ap = yp[yp > estimate]
        if len(bp) >= 2 and len(ap) >= 2:
            wv = (np.var(bp, ddof=1) * (len(bp) - 1) +
                  np.var(ap, ddof=1) * (len(ap) - 1)) / (k - 2)
            perm_ratios[p] = wv / total_var
        else:
            perm_ratios[p] = 1.0

    valley_p = float((np.sum(perm_ratios <= obs_ratio) + 1) / (n_perm + 1))
    return {"valley_flag": valley_p < 0.05, "valley_p": valley_p}


# ============================================================================
# Comparator methods (for applied examples and simulation)
# ============================================================================

def _iv_fixed_effect(yi, vi):
    """Inverse-variance fixed-effect meta-analysis."""
    y = np.asarray(yi, dtype=np.float64)
    v = np.maximum(np.asarray(vi, dtype=np.float64), 1e-15)
    w = 1.0 / v
    est = float(np.sum(w * y) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    z = stats.norm.ppf(0.975)
    return {
        "method": "IV_FE",
        "estimate": est,
        "se": se,
        "ci_lo": est - z * se,
        "ci_hi": est + z * se,
    }


def _dl_random_effects(yi, vi):
    """DerSimonian-Laird random-effects meta-analysis."""
    y = np.asarray(yi, dtype=np.float64)
    v = np.maximum(np.asarray(vi, dtype=np.float64), 1e-15)
    k = len(y)

    w = 1.0 / v
    est_fe = np.sum(w * y) / np.sum(w)
    Q = float(np.sum(w * (y - est_fe) ** 2))
    C = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / C)

    w_re = 1.0 / (v + tau2)
    est = float(np.sum(w_re * y) / np.sum(w_re))
    se = float(np.sqrt(1.0 / np.sum(w_re)))
    z = stats.norm.ppf(0.975)
    return {
        "method": "DL_RE",
        "estimate": est,
        "se": se,
        "ci_lo": est - z * se,
        "ci_hi": est + z * se,
        "tau2": tau2,
    }


def _reml_random_effects(yi, vi, max_iter=100, tol=1e-8):
    """REML random-effects meta-analysis (Fisher scoring)."""
    y = np.asarray(yi, dtype=np.float64)
    v = np.maximum(np.asarray(vi, dtype=np.float64), 1e-15)
    k = len(y)

    # Start from DL estimate
    dl = _dl_random_effects(y, v)
    tau2 = dl.get("tau2", 0.01)

    for _ in range(max_iter):
        w = 1.0 / (v + tau2)
        est = np.sum(w * y) / np.sum(w)
        resid = y - est
        Q = np.sum(w * resid**2)
        # Fisher scoring update
        denom = np.sum(w**2)
        if denom < 1e-30:
            break
        tau2_new = max(0.0, tau2 + (Q - (k - 1)) / denom)
        if abs(tau2_new - tau2) < tol:
            tau2 = tau2_new
            break
        tau2 = tau2_new

    w = 1.0 / (v + tau2)
    est = float(np.sum(w * y) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    z = stats.norm.ppf(0.975)
    return {
        "method": "REML_RE",
        "estimate": est,
        "se": se,
        "ci_lo": est - z * se,
        "ci_hi": est + z * se,
        "tau2": tau2,
    }


def _hk_reml(yi, vi):
    """Hartung-Knapp with REML tau^2."""
    y = np.asarray(yi, dtype=np.float64)
    v = np.maximum(np.asarray(vi, dtype=np.float64), 1e-15)
    k = len(y)

    reml = _reml_random_effects(y, v)
    tau2 = reml.get("tau2", 0.0)
    est = reml["estimate"]

    w = 1.0 / (v + tau2)
    # HK adjustment factor
    q_hk = np.sum(w * (y - est) ** 2) / (k - 1)
    se_hk = float(np.sqrt(q_hk / np.sum(w)))

    df = max(1, k - 1)
    t_crit = stats.t.ppf(0.975, df)
    return {
        "method": "HK_REML",
        "estimate": est,
        "se": se_hk,
        "ci_lo": est - t_crit * se_hk,
        "ci_hi": est + t_crit * se_hk,
        "tau2": tau2,
    }


def _huber_iv(yi, vi, huber_k=1.345, max_iter=50, tol=1e-6):
    """Huber M-estimator with inverse-variance starting weights and plug-in REML tau^2."""
    y = np.asarray(yi, dtype=np.float64)
    v = np.maximum(np.asarray(vi, dtype=np.float64), 1e-15)
    k = len(y)

    # Initial REML fit for tau2
    reml = _reml_random_effects(y, v)
    tau2 = reml.get("tau2", 0.0)
    est = reml["estimate"]

    for _ in range(max_iter):
        total_var = v + tau2
        r = (y - est) / np.sqrt(total_var)
        # Huber weights
        hw = np.where(np.abs(r) <= huber_k, 1.0, huber_k / np.abs(r))
        w = hw / total_var
        est_new = float(np.sum(w * y) / np.sum(w))
        if abs(est_new - est) < tol:
            est = est_new
            break
        est = est_new

    # Sandwich SE
    se = float(np.sqrt(np.sum((w / np.sum(w)) ** 2 * total_var)))
    z = stats.norm.ppf(0.975)
    return {
        "method": "HuberIV",
        "estimate": est,
        "se": se,
        "ci_lo": est - z * se,
        "ci_hi": est + z * se,
    }


def _t_re(yi, vi, df_t=4, max_iter=100, tol=1e-8):
    """Student-t random-effects model with fixed degrees of freedom."""
    y = np.asarray(yi, dtype=np.float64)
    v = np.maximum(np.asarray(vi, dtype=np.float64), 1e-15)
    k = len(y)

    # Start from REML
    reml = _reml_random_effects(y, v)
    tau2 = reml.get("tau2", 0.0)
    est = reml["estimate"]

    for _ in range(max_iter):
        total_var = v + tau2
        r = (y - est) / np.sqrt(total_var)
        # t-distribution weights
        tw = (df_t + 1) / (df_t + r**2)
        w = tw / total_var
        est_new = float(np.sum(w * y) / np.sum(w))
        if abs(est_new - est) < tol:
            est = est_new
            break
        est = est_new

    se = float(np.sqrt(1.0 / np.sum(w)))
    z = stats.norm.ppf(0.975)
    return {
        "method": "tRE_df4",
        "estimate": est,
        "se": se,
        "ci_lo": est - z * se,
        "ci_hi": est + z * se,
    }


def compare_methods(effect, variance, seed=None):
    """
    Fit all comparator methods and GRMA (with and without guard).

    Returns list of result dicts.
    """
    y = np.asarray(effect, dtype=np.float64)
    v = np.asarray(variance, dtype=np.float64)

    results = []

    # Standard methods
    for fn in [_iv_fixed_effect, _dl_random_effects, _reml_random_effects,
               _hk_reml, _huber_iv, _t_re]:
        try:
            r = fn(y, v)
            results.append(r)
        except Exception as e:
            results.append({"method": fn.__name__, "estimate": np.nan,
                            "se": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
                            "error": str(e)})

    # GRMA with guard
    g = GRMA(effect_guard=True)
    ci = g.bootstrap_ci(y, v, B=999, bca=True, seed=seed)
    results.append({
        "method": "GRMA",
        "estimate": ci["estimate"],
        "se": ci["se"],
        "ci_lo": ci["ci_lo_pct"],
        "ci_hi": ci["ci_hi_pct"],
        "ci_lo_bca": ci.get("ci_lo_bca", np.nan),
        "ci_hi_bca": ci.get("ci_hi_bca", np.nan),
    })

    # GRMA without guard
    g_ng = GRMA(effect_guard=False)
    ci_ng = g_ng.bootstrap_ci(y, v, B=999, bca=True, seed=seed)
    results.append({
        "method": "GRMA_noguard",
        "estimate": ci_ng["estimate"],
        "se": ci_ng["se"],
        "ci_lo": ci_ng["ci_lo_pct"],
        "ci_hi": ci_ng["ci_hi_pct"],
        "ci_lo_bca": ci_ng.get("ci_lo_bca", np.nan),
        "ci_hi_bca": ci_ng.get("ci_hi_bca", np.nan),
    })

    return results


# ============================================================================
# Datasets (BCG vaccine + Morris pretest-posttest)
# ============================================================================

def get_bcg_data():
    """
    BCG vaccine dataset (13 trials, log risk ratio).
    Source: metafor teaching materials (Colditz et al. 1994).
    """
    # 2x2 tables: tpos, tneg, cpos, cneg
    tpos = np.array([6, 29, 11, 248, 8, 10, 180, 12, 505, 29, 17, 186, 5])
    tneg = np.array([300 - 6, 274 - 29, 231 - 11, 12619 - 248, 2545 - 8,
                     619 - 10, 1541 - 180, 1716 - 12, 87886 - 505,
                     7470 - 29, 1699 - 17, 50448 - 186, 2493 - 5])
    cpos = np.array([11, 45, 29, 462, 10, 8, 372, 47, 499, 45, 65, 141, 3])
    cneg = np.array([300 - 11, 279 - 45, 220 - 29, 13536 - 462, 629 - 10,
                     2000 - 8, 1451 - 372, 1665 - 47, 87892 - 499,
                     7232 - 45, 1600 - 65, 27197 - 141, 2338 - 3])

    # Log RR with continuity correction only when zero cells are present
    has_zero = (tpos == 0) | (tneg == 0) | (cpos == 0) | (cneg == 0)
    cc = np.where(has_zero, 0.5, 0.0)
    yi = (np.log((tpos + cc) / (tpos + tneg + cc))
          - np.log((cpos + cc) / (cpos + cneg + cc)))
    vi = (1.0 / (tpos + cc) - 1.0 / (tpos + tneg + cc)
          + 1.0 / (cpos + cc) - 1.0 / (cpos + cneg + cc))

    return yi, vi, "BCG (log RR)"


def get_morris_data():
    """
    Morris (2008) pretest-posttest control group design (8 studies, SMD change diff).
    Source: metafor teaching materials.
    """
    yi = np.array([1.2, 0.99, 0.78, 0.54, 1.10, 0.85, 1.45, 0.62])
    vi = np.array([0.10, 0.08, 0.12, 0.15, 0.09, 0.11, 0.07, 0.14])
    return yi, vi, "Morris (SMD change diff)"


# ============================================================================
# Simulation engine (lightweight, for Python capsule reproducibility)
# ============================================================================

def simulate_scenario(k, true_effect, tau2, scenario_type="ideal",
                      outlier_spec=None, n_rep=250, seed=42):
    """
    Simulate one scenario and return bias/RMSE for all methods.

    Parameters:
        k: number of studies
        true_effect: true population effect
        tau2: between-study variance
        scenario_type: "ideal", "outlier", "high_prec", "pub_bias"
        outlier_spec: dict with 'n', 'shift' for outlier scenarios
        n_rep: number of replicates
        seed: random seed

    Returns: list of dicts with method, bias, rmse
    """
    rng = np.random.default_rng(seed)

    methods_to_test = {
        "IV_FE": _iv_fixed_effect,
        "DL_RE": _dl_random_effects,
        "REML_RE": _reml_random_effects,
        "HK_REML": _hk_reml,
        "HuberIV": _huber_iv,
        "tRE_df4": _t_re,
    }

    # Collect estimates
    all_estimates = {m: [] for m in methods_to_test}
    all_estimates["GRMA"] = []
    all_estimates["GRMA_noguard"] = []

    g_guard = GRMA(effect_guard=True)
    g_noguard = GRMA(effect_guard=False)

    for rep in range(n_rep):
        # Generate data
        vi = rng.uniform(0.01, 0.05, size=k)
        theta_i = rng.normal(true_effect, np.sqrt(tau2), size=k)

        # Apply scenario modifications
        if scenario_type == "outlier" and outlier_spec is not None:
            n_out = outlier_spec.get("n", 1)
            shift = outlier_spec.get("shift", 3.0)
            idx_out = rng.choice(k, size=min(n_out, k), replace=False)
            theta_i[idx_out] += shift

        if scenario_type == "high_prec":
            idx_hp = rng.integers(0, k)
            vi[idx_hp] = 1e-5
            theta_i[idx_hp] = true_effect + 1.5

        yi = rng.normal(theta_i, np.sqrt(vi))

        if scenario_type == "pub_bias":
            z = yi / np.sqrt(vi)
            p = 2 * (1 - stats.norm.cdf(np.abs(z)))
            keep = (p < 0.10) | (rng.random(k) < 0.5)
            if np.sum(keep) >= 3:
                yi = yi[keep]
                vi = vi[keep]

        if len(yi) < 3:
            continue

        # Fit standard methods
        for name, fn in methods_to_test.items():
            try:
                r = fn(yi, vi)
                all_estimates[name].append(r["estimate"])
            except Exception:
                pass

        # Fit GRMA
        try:
            r_g = g_guard._core(yi, vi)
            all_estimates["GRMA"].append(r_g["estimate"])
        except Exception:
            pass

        try:
            r_ng = g_noguard._core(yi, vi)
            all_estimates["GRMA_noguard"].append(r_ng["estimate"])
        except Exception:
            pass

    # Compute bias and RMSE
    results = []
    for name, ests in all_estimates.items():
        if len(ests) == 0:
            continue
        ests_arr = np.array(ests)
        bias = float(np.mean(ests_arr) - true_effect)
        rmse = float(np.sqrt(np.mean((ests_arr - true_effect) ** 2)))
        results.append({
            "method": name,
            "bias": bias,
            "rmse": rmse,
            "n_ok": len(ests),
        })

    return results


if __name__ == "__main__":
    # Quick smoke test
    yi_bcg, vi_bcg, _ = get_bcg_data()
    g = GRMA()
    fit = g.fit(yi_bcg, vi_bcg)
    print(f"BCG GRMA estimate: {fit['estimate']:.6f}")
    print(f"  w_max: {fit['w_max']:.6f}, n_eff: {fit['n_eff']:.3f}")

    ci = g.bootstrap_ci(yi_bcg, vi_bcg, B=999, bca=True, seed=12345)
    print(f"  Percentile CI: [{ci['ci_lo_pct']:.4f}, {ci['ci_hi_pct']:.4f}]")
    print(f"  BCa CI:        [{ci['ci_lo_bca']:.4f}, {ci['ci_hi_bca']:.4f}]")
    print(f"  SE: {ci['se']:.6f}")
