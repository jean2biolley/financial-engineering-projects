import numpy as np
import cvxpy as cp
import pandas as pd
from typing import Tuple, List, Optional, Dict
import market_data

# === Configuration ===
# === Configuration ===
R_F = 0.0362

N_TARGET = 30
K = 3
M_MAX = K / N_TARGET
BETA_MAX = 1.0
W_MIN = 0.005

# Turnover Constraints
TURNOVER_TOTAL_MIN = 0.0
TURNOVER_TOTAL_MAX = 100.0
TURNOVER_BASE_COST_FACTOR = 10.0
TURNOVER_VAR_COST_FACTOR = 500

EPS = 1e-6


def get_geo_indices(tickers: List[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    We classify the assets into geographical regions (US, Europe, Asia/Emerging) and return the indices corresponding to each region.
    """
    us_idx = []
    europe_idx = []
    asia_idx = []

    for i, ticker in enumerate(tickers):
        if any(x in ticker for x in
               ['.BR', '.DE', '.PA', '.MI', '.MC', '.AS', '.CO', '.L', '.SW', '.OL', '.ST', '.WA']):
            europe_idx.append(i)
        elif any(x in ticker for x in ['.T', '.KS', '.AX', '.HK', '.KL', '.SZ', '.SS', '.NS', '.SI', '.SR']):
            asia_idx.append(i)
        elif ticker in ['BBVA', 'BABA', 'BBD', 'SONY', 'TM', 'HMC', 'NTDOY', 'MQBKY', 'NOK', 'VALE', 'MT',
                        'INFY', 'UL', 'ACGBY', 'BP', 'FUJHY', 'KAEPY', 'ALIZF']:
            if ticker in ['SONY', 'TM', 'HMC', 'NTDOY', 'MQBKY', 'FUJHY', 'KAEPY', 'INFY', 'BABA']:
                asia_idx.append(i)
            elif ticker in ['BBVA', 'NOK', 'UL', 'ACGBY', 'BP', 'MT', 'ALIZF']:
                europe_idx.append(i)
            elif ticker in ['BBD', 'VALE']:
                asia_idx.append(i)
        else:
            us_idx.append(i)

    return np.array(us_idx), np.array(europe_idx), np.array(asia_idx)


def load_covariance_matrix(file_path: str, sheet_name: str = "cov5") -> np.ndarray:
    """We load the covariance matrix from an Excel file and ensure it is symmetric."""
    cov_df = pd.read_excel(file_path, sheet_name=sheet_name, index_col=0)
    sigma = cov_df.to_numpy()
    sigma = (sigma + sigma.T) / 2
    return sigma


def solve_optimal_portfolio(
        mu: np.ndarray,
        sigma: np.ndarray,
        beta: np.ndarray,
        weights0: np.ndarray,
        geo_indices: Tuple[np.ndarray, np.ndarray, np.ndarray],
        risk_aversion: float = 2.0
) -> Optional[np.ndarray]:
    """
    We solve the Mean-Variance Optimization problem using Mixed-Integer Quadratic Programming (MIQP).
    We maximize the quadratic utility (Return - Gamma * Risk) subject to constraints on allocation, budget, and turnover costs.
    Then, the code returns the optimal weights if a solution is found, otherwise None.
    """
    n = len(mu)
    us_idx, europe_idx, asia_idx = geo_indices

    # Decision Variables
    w = cp.Variable(n, nonneg=True)
    is_changed = cp.Variable(n, boolean=True)
    diff = cp.Variable(n, nonneg=True)

    # Constants for turnover
    big_M = 1.0  # Weights are <= 1
    min_trade = 1e-6

    constraints = [
        # Budget and Allocation
        cp.sum(w) == 1,
        w <= M_MAX,
        beta @ w <= BETA_MAX,
        # We apply Geographical Constraints
        cp.sum(w[us_idx]) >= 0.35,
        cp.sum(w[us_idx]) <= 0.45,
        cp.sum(w[europe_idx]) >= 0.25,
        cp.sum(w[europe_idx]) <= 0.35,
        cp.sum(w[asia_idx]) >= 0.25,
        cp.sum(w[asia_idx]) <= 0.35,

        # Turnover Constraints (MIP)
        # |w - w0| <= diff
        w - weights0 <= diff,
        weights0 - w <= diff,

        # diff <= M * indicator
        diff <= big_M * is_changed,
        # diff >= min_trade * indicator (to force indicator=1 if diff>0)
        diff >= min_trade * is_changed,

        # Total Cost Constraint
        # Cost = Base * Sum(Indicator) + Var * Sum(diff) <= 100
        TURNOVER_BASE_COST_FACTOR * cp.sum(is_changed) + TURNOVER_VAR_COST_FACTOR * cp.sum(diff) <= TURNOVER_TOTAL_MAX
    ]

    # Objective: Maximize Utility (Return - Risk_Aversion * Sample_Variance)
    # Using quad_form for variance: w.T @ Sigma @ w
    risk_term = cp.quad_form(w, sigma)
    return_term = mu @ w

    objective = cp.Maximize(return_term - risk_aversion * risk_term)

    problem = cp.Problem(objective, constraints)

    # Solve with available MIQP solvers
    solved = False
    for solver in [cp.CPLEX, cp.GUROBI, cp.CBC, cp.GLPK_MI, cp.SCIPY]:
        if solver in cp.installed_solvers():
            try:
                # SCIPY generally doesn't support MIQP well, but we try as fallback
                problem.solve(solver=solver)
                if problem.status in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
                    solved = True
                    break
            except Exception:
                continue

    if not solved or w.value is None:
        return None

    w_final = np.maximum(w.value, 0)
    w_final /= np.sum(w_final)  # Normalize to ensure sum is exactly 1
    return w_final


def maximize_sharpe_with_turnover(
        mu: np.ndarray,
        sigma: np.ndarray,
        beta: np.ndarray,
        weights0: np.ndarray,
        geo_indices: Tuple[np.ndarray, np.ndarray, np.ndarray]
) -> np.ndarray:
    """
    We perform a grid search over the risk aversion parameter (Gamma) to find the portfolio
    that maximizes the Sharpe Ratio while ensuring all constraints, including turnover, are satisfied.
    """
    print("Scanning Efficient Frontier to Maximize Sharpe with Turnover <= 100...")

    # Log-spaced gammas from very low (risk-seeking/max return) to high (min vol)
    # We want to cover the tangent point region.
    gammas = np.logspace(-1, 2, 20)  # 0.1 to 100

    best_sharpe = -np.inf
    best_w = weights0.copy()

    # Also evaluate initial weights
    risk0 = np.sqrt(weights0.T @ sigma @ weights0)
    ret0 = mu @ weights0
    sharpe0 = (ret0 - R_F) / risk0 if risk0 > 0 else -np.inf
    print(f"  > Initial Portfolio: Sharpe = {sharpe0:.4f}")
    if sharpe0 > best_sharpe:
        best_sharpe = sharpe0
        best_w = weights0

    for gamma in gammas:
        w_res = solve_optimal_portfolio(mu, sigma, beta, weights0, geo_indices, risk_aversion=gamma)

        if w_res is not None:
            # Calculate metrics
            risk = np.sqrt(w_res.T @ sigma @ w_res)
            ret = mu @ w_res
            sharpe = (ret - R_F) / risk if risk > 0 else -np.inf
            # print(f"    Gamma {gamma:.2f}: Sharpe={sharpe:.4f}")

            if sharpe > best_sharpe:
                best_sharpe = sharpe
                best_w = w_res

    print(f"  > Best Found Sharpe: {best_sharpe:.4f}")
    return best_w


def turnover_cost_np(weights: np.ndarray, weights0: np.ndarray, tol: float = 1e-6) -> Tuple[float, float, float]:
    """We calculate the actual turnover cost components (Base + Variable) for reporting purposes."""
    diff = np.abs(weights - weights0)
    changed_mask = diff > tol
    base_cost = TURNOVER_BASE_COST_FACTOR * np.sum(changed_mask)
    variable_cost = np.sum(diff * TURNOVER_VAR_COST_FACTOR)
    return base_cost + variable_cost, base_cost, variable_cost


def print_report(
        w_final: np.ndarray,
        mu: np.ndarray,
        sigma: np.ndarray,
        beta: np.ndarray,
        geo_indices: Tuple[np.ndarray, np.ndarray, np.ndarray],
        turnover_stats: Tuple[float, float, float],
        tickers: List[str]
):
    """
    We generate and print a detailed report of the final portfolio.
    This includes geographical allocation, market cap allocation, risk/return metrics, and top holdings.
    """
    us_idx, europe_idx, asia_idx = geo_indices
    total_cost, base_cost, var_cost = turnover_stats

    ret = float(mu @ w_final)
    risk = float(np.sqrt(w_final.T @ sigma @ w_final))
    sharpe = (ret - R_F) / risk if risk > 0 else np.nan

    n_selected = np.sum(w_final > EPS)

    print("\n=== Final Portfolio Report ===")
    print(f"Geographical Allocation:")
    print(f"  US:     {np.sum(w_final[us_idx]) * 100:.2f}%")
    print(f"  Europe: {np.sum(w_final[europe_idx]) * 100:.2f}%")
    print(f"  Asia:   {np.sum(w_final[asia_idx]) * 100:.2f}%")

    print(f"\nMetrics:")
    print(f"  Result Beta: {np.dot(beta, w_final):.4f} (Max: {BETA_MAX})")
    print(f"  Max Weight: {np.max(w_final) * 100:.2f}% (Limit: {M_MAX * 100:.2f}%)")
    print(f"  Assets Selected: {n_selected}")
    print(f"  Turnover Cost: {total_cost:.2f} (Base: {base_cost:.2f}, Var: {var_cost:.2f})")
    print(f"  Exp Return: {ret * 100:.4f}% | Risk: {risk * 100:.4f}% | Sharpe: {sharpe:.4f}")

    print(f"\nTop Holdings:")
    for i in np.argsort(w_final)[-15:][::-1]:
        if w_final[i] > EPS:
            print(f"  {tickers[i]:10s}: {w_final[i] * 100:7.2f}%")


def export_results(w_final: np.ndarray, mu: np.ndarray, beta: np.ndarray, tickers: List[str],
                   geo_indices: Tuple[np.ndarray, np.ndarray, np.ndarray]):
    """
    We export the optimization results to an Excel file.
    We create two sheets: one with all assets and one with only the selected assets.
    """
    print(f"\n=== Exporting to Excel ===")
    us_idx, europe_idx, _ = geo_indices

    regions = []
    for i in range(len(tickers)):
        if i in us_idx:
            regions.append('US')
        elif i in europe_idx:
            regions.append('Europe')
        else:
            regions.append('Asia')

    df = pd.DataFrame({
        'Ticker': tickers,
        'Weight': w_final,
        'Weight_Percent': w_final * 100,
        'Region': regions,

        'Beta': beta,
        'Exp_Return': mu,
        'Selected': (w_final > EPS).astype(int)
    })

    selected_df = df[df['Selected'] == 1].drop(columns=['Selected'])

    filename = 'portfolio_weights5.xlsx'
    try:
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            df.to_excel(writer, index=False, sheet_name='All Assets')
            selected_df.to_excel(writer, index=False, sheet_name='Selected Assets Only')
        print(f"Exported to {filename}")
    except Exception as e:
        print(f"Export failed: {e}")


def main():
    # 1. Load Data
    print("Loading data...")
    market_data_dict = market_data.load_market_data("main techno.xlsx")

    mu = market_data_dict["mu"]
    beta = market_data_dict["beta"]

    weights0 = market_data_dict["weights0"]
    tickers = market_data_dict["tickers"]

    # 2. We set up the indices for regions and market caps, and load the covariance matrix.
    geo_indices = get_geo_indices(tickers)

    # Check dimensions
    if len(mu) != len(tickers):
        raise ValueError(f"Mismatch: {len(tickers)} tickers but {len(mu)} returns.")

    sigma = load_covariance_matrix("covperso.xlsx")

    # Ensure Sigma matches dimensions (if not, slice it)
    if sigma.shape[0] != len(tickers):
        print(f"Warning: Sigma shape {sigma.shape} does not match assets {len(tickers)}. Slicing Sigma.")
        sigma = sigma[:len(tickers), :len(tickers)]

    print("\n=== Starting Optimization ===")
    # 3. Optimize Sharpe with Turnover Constraints (Grid Search)
    w_final = maximize_sharpe_with_turnover(mu, sigma, beta, weights0, geo_indices)

    # Check final constraints
    turnover_stats = turnover_cost_np(w_final, weights0)
    if (turnover_stats[0] > TURNOVER_TOTAL_MAX + 1e-5 or
            turnover_stats[0] < TURNOVER_TOTAL_MIN - 1e-5):
        print("WARNING: Final portfolio violates turnover constraints.")

    # 6. We report the performance metrics and export the final portfolio to Excel.
    print_report(w_final, mu, sigma, beta, geo_indices, turnover_stats, tickers)
    export_results(w_final, mu, beta, tickers, geo_indices)


if __name__ == "__main__":
    main()