"""
Black-Scholes option pricer with Greeks, a put-call parity check and a Monte Carlo check.

Run `python my_pricer.py` to see all the checks.
"""
import numpy as np
from scipy.stats import norm


def d1_d2(S, K, T, r, sigma):
    """The two quantities that appear in every Black-Scholes formula."""
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return d1, d2


def bs_price(S, K, T, r, sigma, kind="call"):
    """Black-Scholes price of a European call or put (no dividends)."""
    d1, d2 = d1_d2(S, K, T, r, sigma)
    if kind == "call":
        return S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    return K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)


def bs_delta_gamma_vega(S, K, T, r, sigma):
    """Call delta, gamma and vega (vega is per 1.00 change in volatility)."""
    d1, _ = d1_d2(S, K, T, r, sigma)
    delta = norm.cdf(d1)
    gamma = norm.pdf(d1) / (S * sigma * np.sqrt(T))
    vega = S * norm.pdf(d1) * np.sqrt(T)
    return delta, gamma, vega


def mc_call(S, K, T, r, sigma, n=200_000, seed=0):
    """Monte Carlo call price. Returns (price, standard_error)."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n)
    ST = S * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * z)
    payoff = np.exp(-r * T) * np.maximum(ST - K, 0.0)
    return payoff.mean(), payoff.std(ddof=1) / np.sqrt(n)


if __name__ == "__main__":
    S, K, T, r, sigma = 100, 100, 1, 0.05, 0.2

    call = bs_price(S, K, T, r, sigma, "call")
    put = bs_price(S, K, T, r, sigma, "put")
    print(f"Call: {call:.4f}   Put: {put:.4f}")

    delta, gamma, vega = bs_delta_gamma_vega(S, K, T, r, sigma)
    print(f"Delta: {delta:.4f}   Gamma: {gamma:.4f}   Vega: {vega:.4f}")

    # Put-call parity: C - P should equal S - K*exp(-rT)
    print(f"Parity check: C - P = {call - put:.6f}, S - K*exp(-rT) = {S - K * np.exp(-r * T):.6f}")

    # Monte Carlo should land within a few standard errors of Black-Scholes
    mc, se = mc_call(S, K, T, r, sigma)
    print(f"Monte Carlo call: {mc:.4f} (+/- {se:.4f}), Black-Scholes: {call:.4f}")
