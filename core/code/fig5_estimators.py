from __future__ import annotations

import argparse

import hashlib

import importlib.util

import itertools

import json

import math

import platform

import sys

from pathlib import Path

import numpy as np

import pandas as pd

import scipy

import sklearn

from joblib import Parallel, delayed

from scipy.stats import norm, t as student_t

from sklearn.ensemble import RandomForestRegressor

from sklearn.metrics import mean_squared_error, r2_score

from sklearn.model_selection import GroupKFold

COHORTS = {
    "optionC_date_aligned_n107_order_locked": 107,
    "optionC_observed_RGI_n75_order_locked": 75,
}

GCM_GROUP_SEED = 1

GCOMP_GROUP_SEED = 7

KNOCKOFF_GROUP_SEED = 17

N_SPLITS = 5

KNOCKOFF_SEED_COUNT = 30

KNOCKOFF_TREES = 600

GCM_TREES = 200

GCOMP_TREES = 500

GCOMP_BOOTSTRAPS = 2000

def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result

def grouped_splits(
    frame: pd.DataFrame,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    groups = frame["site"].astype(str).to_numpy()
    splitter = GroupKFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=seed,
    )
    splits = list(splitter.split(frame, groups=groups))
    for train, test in splits:
        overlap = set(groups[train]) & set(groups[test])
        if overlap:
            raise AssertionError(f"Site leakage in grouped split: {overlap}")
    return splits

def crossfit_residuals_grouped(
    predictors: np.ndarray,
    outcome: np.ndarray,
    splits: list[tuple[np.ndarray, np.ndarray]],
    n_jobs: int = 4,
) -> np.ndarray:
    if predictors.shape[1] == 0:
        return outcome - outcome.mean()
    residuals = np.zeros(len(outcome), dtype=float)
    assigned = np.zeros(len(outcome), dtype=int)
    for train, test in splits:
        model = RandomForestRegressor(
            GCM_TREES,
            min_samples_leaf=3,
            random_state=0,
            n_jobs=n_jobs,
        ).fit(predictors[train], outcome[train])
        residuals[test] = outcome[test] - model.predict(predictors[test])
        assigned[test] += 1
    if not np.all(assigned == 1):
        raise AssertionError("Grouped GCM residuals were not assigned exactly once")
    return residuals

def cluster_statistics(
    products: np.ndarray,
    sites: np.ndarray,
) -> dict[str, float]:
    n = len(products)
    mean_product = float(products.mean())
    site_frame = pd.DataFrame(
        {"site": sites.astype(str), "product": products}
    )
    site_sums = site_frame.groupby("site", sort=True)["product"].sum()
    site_counts = site_frame.groupby("site", sort=True)["product"].size()
    clusters = len(site_sums)
    centered_sums = site_sums.to_numpy() - site_counts.to_numpy() * mean_product
    se_cluster = (
        math.sqrt(clusters / (clusters - 1) * np.sum(centered_sums**2)) / n
    )
    t_cluster = mean_product / se_cluster if se_cluster > 0 else np.nan
    p_cluster = (
        2.0 * student_t.sf(abs(t_cluster), df=clusters - 1)
        if np.isfinite(t_cluster)
        else np.nan
    )
    signs = np.asarray(list(itertools.product([-1.0, 1.0], repeat=clusters)))
    sign_flipped = signs @ site_sums.to_numpy(dtype=float)
    all_plus_index = int(
        np.flatnonzero((signs == 1.0).all(axis=1))[0]
    )
    observed_reference = float(abs(sign_flipped[all_plus_index]))
    comparison_floor = np.nextafter(observed_reference, -np.inf)
    p_signflip = float(
        np.mean(np.abs(sign_flipped) >= comparison_floor)
    )
    attainable_floor = 2.0 / (2**clusters)
    if p_signflip < attainable_floor or p_signflip == 0.0:
        raise AssertionError("Invalid exact site sign-flip P below attainable floor")
    return {
        "mean_residual_product": mean_product,
        "n_site_clusters": clusters,
        "cluster_df": clusters - 1,
        "se_cluster": se_cluster,
        "T_cluster": t_cluster,
        "p_t_cluster_df": p_cluster,
        "p_exact_site_signflip_2s": p_signflip,
        "site_signflip_enumerations": 2**clusters,
        "site_signflip_observed_reference_abs_sum": observed_reference,
    }

def gcm_driver_worker(
    cohort: str,
    driver_spec: tuple[str, str, str],
    subset: pd.DataFrame,
    community: np.ndarray,
    outcome_residuals: dict[int, np.ndarray],
    splits: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Compute all 11 conditioning depths for one GCM driver."""
    label, column, group = driver_spec
    exposure = subset[column].to_numpy(dtype=float)
    sites = subset["site"].astype(str).to_numpy()
    result_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    for k in range(11):
        exposure_residuals = crossfit_residuals_grouped(
            community[:, :k],
            exposure,
            splits,
            n_jobs=1,
        )
        products = exposure_residuals * outcome_residuals[k]
        variance = float((products**2).mean() - products.mean() ** 2)
        t_normal = (
            math.sqrt(len(subset)) * float(products.mean()) / math.sqrt(variance)
            if variance > 0
            else np.nan
        )
        row = {
            "cohort": cohort,
            "driver": label,
            "column": column,
            "group": group,
            "n": len(subset),
            "community_pcs_k": k,
            "crossfit": "5-fold GroupKFold, sites held out",
            "group_split_seed": GCM_GROUP_SEED,
            "rf_trees": GCM_TREES,
            "rf_min_samples_leaf": 3,
            "rf_random_state": 0,
            "T_normal_comparability": t_normal,
            "p_normal_2s_comparability": (
                2.0 * (1.0 - norm.cdf(abs(t_normal)))
                if np.isfinite(t_normal)
                else np.nan
            ),
        }
        row.update(cluster_statistics(products, sites))
        result_rows.append(row)
        product_frame = pd.DataFrame(
            {"site": sites, "residual_product": products}
        )
        for site, site_frame in product_frame.groupby("site", sort=True):
            score_rows.append(
                {
                    "cohort": cohort,
                    "driver": label,
                    "column": column,
                    "community_pcs_k": k,
                    "site": site,
                    "site_n": len(site_frame),
                    "site_residual_product_sum": (
                        site_frame["residual_product"].sum()
                    ),
                    "site_residual_product_mean": (
                        site_frame["residual_product"].mean()
                    ),
                }
            )
    return result_rows, score_rows

def run_gcm_site_blocked(
    source,
    cohorts: dict[str, pd.DataFrame],
) -> None:
    rows: list[dict[str, object]] = []
    all_score_rows: list[dict[str, object]] = []
    required = [
        column for _, column, _ in source.GCM_DRIVERS
    ] + source.PCS + ["ARG_total", "site"]
    for cohort, data in cohorts.items():
        subset = data.dropna(subset=required).reset_index(drop=True)
        splits = grouped_splits(subset, GCM_GROUP_SEED)
        outcome = subset["ARG_total"].to_numpy(dtype=float)
        sites = subset["site"].astype(str).to_numpy()
        community = subset[source.PCS].to_numpy(dtype=float)
        outcome_residuals = {
            k: crossfit_residuals_grouped(community[:, :k], outcome, splits)
            for k in range(11)
        }
        completed = Parallel(n_jobs=5, backend="loky", verbose=5)(
            delayed(gcm_driver_worker)(
                cohort,
                driver_spec,
                subset,
                community,
                outcome_residuals,
                splits,
            )
            for driver_spec in source.GCM_DRIVERS
        )
        for driver_rows, score_rows in completed:
            rows.extend(driver_rows)
            all_score_rows.extend(score_rows)
        print(f"[SITE-BLOCKED GCM] {cohort}: complete", flush=True)
    result = pd.DataFrame(rows)
    for cohort in COHORTS:
        for k in range(11):
            mask = result["cohort"].eq(cohort) & result["community_pcs_k"].eq(k)
            result.loc[mask, "q_bh23_normal_comparability"] = bh_adjust(
                result.loc[mask, "p_normal_2s_comparability"].to_numpy()
            )
            result.loc[mask, "q_bh23_t_cluster_df"] = bh_adjust(
                result.loc[mask, "p_t_cluster_df"].to_numpy()
            )
            result.loc[mask, "q_bh23_site_signflip"] = bh_adjust(
                result.loc[mask, "p_exact_site_signflip_2s"].to_numpy()
            )
    result.to_csv(OUTPUT / "gcm_site_blocked.csv", index=False)
    pd.DataFrame(all_score_rows).to_csv(
        OUTPUT / "gcm_site_cluster_scores.csv",
        index=False,
    )

def cluster_bootstrap_means(
    values: np.ndarray,
    sites: np.ndarray,
    seed: int,
    iterations: int,
) -> np.ndarray:
    unique_sites = np.sort(np.unique(sites.astype(str)))
    by_site = {
        site: values[sites.astype(str) == site]
        for site in unique_sites
    }
    random = np.random.RandomState(seed)
    means = np.empty(iterations, dtype=float)
    for index in range(iterations):
        sampled_sites = random.choice(unique_sites, size=len(unique_sites), replace=True)
        sampled_values = np.concatenate([by_site[site] for site in sampled_sites])
        means[index] = sampled_values.mean()
    return means

def empirical_cdf(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    return np.array([(values <= point).mean() for point in grid])

def run_gcomp_site_blocked(
    source,
    cohorts: dict[str, pd.DataFrame],
) -> None:
    rows: list[dict[str, object]] = []
    for cohort, data in cohorts.items():
        required = source.GCOMP_FEATURES + ["ARG_total", "site"]
        subset = data.dropna(subset=required).reset_index(drop=True)
        predictors = subset[source.GCOMP_FEATURES].copy()
        outcome = subset["ARG_total"].to_numpy(dtype=float)
        sites = subset["site"].astype(str).to_numpy()
        outcome_sd = outcome.std(ddof=0)
        splits = grouped_splits(subset, GCOMP_GROUP_SEED)
        for driver_index, (label, column, mode) in enumerate(
            source.GCOMP_DRIVERS
        ):
            high_predictions = np.zeros(len(subset), dtype=float)
            low_predictions = np.zeros(len(subset), dtype=float)
            base_predictions = np.zeros(len(subset), dtype=float)
            assigned = np.zeros(len(subset), dtype=int)
            global_low = (
                float(np.quantile(predictors[column], 0.1))
                if mode == "quantile"
                else 0.0
            )
            global_high = (
                float(np.quantile(predictors[column], 0.9))
                if mode == "quantile"
                else 1.0
            )
            for train, test in splits:
                model = RandomForestRegressor(
                    GCOMP_TREES,
                    min_samples_leaf=3,
                    random_state=0,
                    n_jobs=4,
                ).fit(predictors.iloc[train], outcome[train])
                high = predictors.iloc[test].copy()
                low = predictors.iloc[test].copy()
                high[column] = global_high
                low[column] = global_low
                base_predictions[test] = model.predict(predictors.iloc[test])
                high_predictions[test] = model.predict(high)
                low_predictions[test] = model.predict(low)
                assigned[test] += 1
            if not np.all(assigned == 1):
                raise AssertionError("Grouped G-computation predictions incomplete")
            individual_shift = high_predictions - low_predictions
            random = np.random.RandomState(700 + driver_index)
            row_bootstrap = np.empty(GCOMP_BOOTSTRAPS, dtype=float)
            for bootstrap_index in range(GCOMP_BOOTSTRAPS):
                indices = random.randint(0, len(subset), len(subset))
                row_bootstrap[bootstrap_index] = individual_shift[indices].mean()
            site_bootstrap = cluster_bootstrap_means(
                individual_shift,
                sites,
                seed=1700 + driver_index,
                iterations=GCOMP_BOOTSTRAPS,
            )
            grid = np.linspace(
                min(high_predictions.min(), low_predictions.min()),
                max(high_predictions.max(), low_predictions.max()),
                400,
            )
            delta_cdf = np.max(
                np.abs(
                    empirical_cdf(high_predictions, grid)
                    - empirical_cdf(low_predictions, grid)
                )
            )
            mean_shift = float(individual_shift.mean())
            row_low, row_high = np.percentile(row_bootstrap, [2.5, 97.5])
            site_low, site_high = np.percentile(site_bootstrap, [2.5, 97.5])
            rows.append(
                {
                    "cohort": cohort,
                    "driver": label,
                    "column": column,
                    "mode": mode,
                    "n": len(subset),
                    "n_site_clusters": len(np.unique(sites)),
                    "low_value": global_low,
                    "high_value": global_high,
                    "crossfit": "5-fold GroupKFold, sites held out",
                    "group_split_seed": GCOMP_GROUP_SEED,
                    "rf_trees": GCOMP_TREES,
                    "rf_min_samples_leaf": 3,
                    "rf_random_state": 0,
                    "mean_shift": mean_shift,
                    "mean_shift_sd": mean_shift / outcome_sd,
                    "row_bootstrap_seed": 700 + driver_index,
                    "row_bootstrap_iterations": GCOMP_BOOTSTRAPS,
                    "mean_shift_ci_low_sd_row_boot": row_low / outcome_sd,
                    "mean_shift_ci_high_sd_row_boot": row_high / outcome_sd,
                    "site_bootstrap_seed": 1700 + driver_index,
                    "site_bootstrap_iterations": GCOMP_BOOTSTRAPS,
                    "mean_shift_ci_low_sd_site_boot": site_low / outcome_sd,
                    "mean_shift_ci_high_sd_site_boot": site_high / outcome_sd,
                    "max_abs_delta_cdf": delta_cdf,
                    "site_heldout_base_rmse": math.sqrt(
                        mean_squared_error(outcome, base_predictions)
                    ),
                    "site_heldout_base_r2": r2_score(outcome, base_predictions),
                }
            )
        print(f"[SITE-BLOCKED G-COMP] {cohort}: complete", flush=True)
    pd.DataFrame(rows).to_csv(
        OUTPUT / "gcomputation_site_blocked.csv",
        index=False,
    )

def gaussian_knockoff_parameters(
    standardized_train: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, float, bool]:
    covariance = np.corrcoef(standardized_train, rowvar=False)
    shrunk = bool(np.linalg.eigvalsh(covariance).min() < 1e-3)
    if shrunk:
        covariance = 0.9 * covariance + 0.1 * np.eye(covariance.shape[0])
    minimum_eigenvalue = float(np.linalg.eigvalsh(covariance).min())
    p_features = standardized_train.shape[1]
    s_value = min(1.0, 2.0 * minimum_eigenvalue)
    s_matrix = np.diag(np.repeat(s_value, p_features))
    covariance_inverse = np.linalg.inv(covariance)
    knockoff_covariance = (
        2.0 * s_matrix - s_matrix @ covariance_inverse @ s_matrix
    )
    knockoff_covariance = (
        knockoff_covariance + knockoff_covariance.T
    ) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(knockoff_covariance)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    covariance_root = eigenvectors @ np.diag(np.sqrt(eigenvalues))
    return (
        covariance_inverse,
        s_matrix,
        minimum_eigenvalue,
        s_value,
        shrunk,
        covariance_root,
    )

def generate_knockoff(
    standardized: np.ndarray,
    covariance_inverse: np.ndarray,
    s_matrix: np.ndarray,
    covariance_root: np.ndarray,
    random: np.random.RandomState,
) -> np.ndarray:
    mean = standardized - standardized @ covariance_inverse @ s_matrix
    return (
        mean
        + random.standard_normal(standardized.shape) @ covariance_root.T
    )

def descriptive_threshold(statistics: np.ndarray) -> float:
    candidates = np.sort(np.abs(statistics[statistics != 0]))
    for candidate in candidates:
        false_count = int((statistics <= -candidate).sum())
        selected_count = max(1, int((statistics >= candidate).sum()))
        if false_count / selected_count <= 0.2:
            return float(candidate)
    return math.inf

def knockoff_fold_worker(
    cohort: str,
    fold_index: int,
    train: np.ndarray,
    test: np.ndarray,
    matrix: np.ndarray,
    outcome: np.ndarray,
    site_values: np.ndarray,
    columns: list[str],
) -> dict[str, object]:
    """Run one outer site-held-out knockoff fold in an isolated process."""
    p_features = len(columns)
    train_mean = matrix[train].mean(axis=0)
    train_sd = matrix[train].std(axis=0, ddof=0)
    if np.any(train_sd == 0):
        zero_columns = np.asarray(columns)[train_sd == 0].tolist()
        raise ValueError(
            f"{cohort} fold {fold_index + 1}: zero-SD columns {zero_columns}"
        )
    train_x = (matrix[train] - train_mean) / train_sd
    test_x = (matrix[test] - train_mean) / train_sd
    (
        covariance_inverse,
        s_matrix,
        minimum_eigenvalue,
        s_value,
        shrunk,
        covariance_root,
    ) = gaussian_knockoff_parameters(train_x)
    w_values = np.zeros((KNOCKOFF_SEED_COUNT, p_features))
    original_delta = np.zeros_like(w_values)
    knockoff_delta = np.zeros_like(w_values)
    performance_rows: list[dict[str, object]] = []
    for seed_index in range(KNOCKOFF_SEED_COUNT):
        random = np.random.RandomState(300 + seed_index)
        train_knockoff = generate_knockoff(
            train_x,
            covariance_inverse,
            s_matrix,
            covariance_root,
            random,
        )
        test_knockoff = generate_knockoff(
            test_x,
            covariance_inverse,
            s_matrix,
            covariance_root,
            random,
        )
        augmented_train = np.column_stack([train_x, train_knockoff])
        augmented_test = np.column_stack([test_x, test_knockoff])
        model = RandomForestRegressor(
            KNOCKOFF_TREES,
            min_samples_leaf=3,
            random_state=seed_index,
            n_jobs=1,
        ).fit(augmented_train, outcome[train])
        baseline = model.predict(augmented_test)
        baseline_mse = mean_squared_error(outcome[test], baseline)
        performance_rows.append(
            {
                "cohort": cohort,
                "fold_one_based": fold_index + 1,
                "seed_index": seed_index,
                "knockoff_seed": 300 + seed_index,
                "rf_random_state": seed_index,
                "train_n": len(train),
                "test_n": len(test),
                "test_sites": ";".join(sorted(np.unique(site_values[test]))),
                "baseline_mse": baseline_mse,
                "baseline_rmse": math.sqrt(baseline_mse),
                "covariance_shrunk": shrunk,
                "minimum_eigenvalue_after_shrinkage": minimum_eigenvalue,
                "equicorrelated_s": s_value,
            }
        )
        for feature_index in range(p_features):
            permutation = np.random.RandomState(
                900000
                + 10000 * fold_index
                + 100 * seed_index
                + feature_index
            ).permutation(len(test))
            permuted_original = augmented_test.copy()
            permuted_original[:, feature_index] = (
                permuted_original[permutation, feature_index]
            )
            permuted_knockoff = augmented_test.copy()
            knockoff_index = p_features + feature_index
            permuted_knockoff[:, knockoff_index] = (
                permuted_knockoff[permutation, knockoff_index]
            )
            original_loss = mean_squared_error(
                outcome[test],
                model.predict(permuted_original),
            )
            knockoff_loss = mean_squared_error(
                outcome[test],
                model.predict(permuted_knockoff),
            )
            original_delta[seed_index, feature_index] = (
                original_loss - baseline_mse
            )
            knockoff_delta[seed_index, feature_index] = (
                knockoff_loss - baseline_mse
            )
            w_values[seed_index, feature_index] = (
                original_loss - knockoff_loss
            )
    return {
        "fold_index": fold_index,
        "test_n": len(test),
        "w_values": w_values,
        "original_delta": original_delta,
        "knockoff_delta": knockoff_delta,
        "minimum_eigenvalue": minimum_eigenvalue,
        "s_value": s_value,
        "shrunk": shrunk,
        "performance_rows": performance_rows,
    }

def run_knockoff_site_blocked(
    source,
    cohorts: dict[str, pd.DataFrame],
) -> None:
    result_frames: list[pd.DataFrame] = []
    performance_rows: list[dict[str, object]] = []
    labels = [item[0] for item in source.KNOCKOFF_ENV] + [
        f"PC{i}" for i in range(1, 11)
    ]
    groups = [item[2] for item in source.KNOCKOFF_ENV] + ["community"] * 10
    columns = [item[1] for item in source.KNOCKOFF_ENV] + source.PCS
    p_features = len(columns)
    for cohort, data in cohorts.items():
        subset = data[["sample_id", "site", "ARG_total"] + columns].dropna()
        subset = subset.reset_index(drop=True)
        matrix = subset[columns].to_numpy(dtype=float)
        outcome = subset["ARG_total"].to_numpy(dtype=float)
        splits = grouped_splits(subset, KNOCKOFF_GROUP_SEED)
        completed = Parallel(n_jobs=N_SPLITS, backend="loky", verbose=5)(
            delayed(knockoff_fold_worker)(
                cohort,
                fold_index,
                train,
                test,
                matrix,
                outcome,
                subset["site"].astype(str).to_numpy(),
                columns,
            )
            for fold_index, (train, test) in enumerate(splits)
        )
        completed = sorted(completed, key=lambda item: item["fold_index"])
        w_values = np.stack([item["w_values"] for item in completed])
        original_delta = np.stack(
            [item["original_delta"] for item in completed]
        )
        knockoff_delta = np.stack(
            [item["knockoff_delta"] for item in completed]
        )
        fold_weights = np.asarray(
            [item["test_n"] for item in completed],
            dtype=float,
        )
        covariance_minima = [
            item["minimum_eigenvalue"] for item in completed
        ]
        s_values = [item["s_value"] for item in completed]
        shrink_count = sum(bool(item["shrunk"]) for item in completed)
        for item in completed:
            performance_rows.extend(item["performance_rows"])
        print(
            f"[SITE-BLOCKED KNOCKOFF] {cohort}: "
            f"{N_SPLITS} folds complete in parallel",
            flush=True,
        )
        normalized_weights = fold_weights / fold_weights.sum()
        mean_w_by_seed = np.tensordot(
            normalized_weights,
            w_values,
            axes=(0, 0),
        )
        mean_original_by_seed = np.tensordot(
            normalized_weights,
            original_delta,
            axes=(0, 0),
        )
        mean_knockoff_by_seed = np.tensordot(
            normalized_weights,
            knockoff_delta,
            axes=(0, 0),
        )
        mean_w = mean_w_by_seed.mean(axis=0)
        threshold = descriptive_threshold(mean_w)
        selected = mean_w >= threshold
        result_frames.append(
            pd.DataFrame(
                {
                    "cohort": cohort,
                    "feature": columns,
                    "label": labels,
                    "group": groups,
                    "n": len(subset),
                    "n_site_clusters": subset["site"].nunique(),
                    "outer_group_split_seed": KNOCKOFF_GROUP_SEED,
                    "knockoff_seed_first": 300,
                    "knockoff_seed_last": 329,
                    "seed_count": KNOCKOFF_SEED_COUNT,
                    "rf_trees": KNOCKOFF_TREES,
                    "rf_min_samples_leaf": 3,
                    "heldout_original_permutation_delta_mse_mean": (
                        mean_original_by_seed.mean(axis=0)
                    ),
                    "heldout_knockoff_permutation_delta_mse_mean": (
                        mean_knockoff_by_seed.mean(axis=0)
                    ),
                    "heldout_W_mean": mean_w,
                    "heldout_W_positive_stability": (
                        mean_w_by_seed > 0
                    ).mean(axis=0),
                    "descriptive_ratio_threshold_0_2": threshold,
                    "selected_descriptive_ratio_0_2": selected,
                    "formal_modelX_fdr_valid": False,
                    "fold_covariance_shrink_count": shrink_count,
                    "fold_minimum_eigenvalue_min": min(covariance_minima),
                    "fold_equicorrelated_s_min": min(s_values),
                }
            )
        )
    pd.concat(result_frames, ignore_index=True).to_csv(
        OUTPUT / "knockoff_site_blocked_predictive_mirror.csv",
        index=False,
    )
    pd.DataFrame(performance_rows).to_csv(
        OUTPUT / "knockoff_site_blocked_fold_performance.csv",
        index=False,
    )
