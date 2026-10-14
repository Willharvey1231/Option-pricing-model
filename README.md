
# Black-Scholes Option Pricer

A small Python implementation of the Black-Scholes model, written to check my understanding of
the maths from my stochastic processes coursework and options trading.

## What it does
- Prices European calls and puts with the Black-Scholes formula
- Computes call delta, gamma and vega
- Checks put-call parity numerically
- Prices a call by Monte Carlo simulation and compares it with the closed-form price


## Example output
For S = 100, K = 100, T = 1 year, r = 5%, sigma = 20%:
```
Call: 10.4506   Put: 5.5735
Delta: 0.6368   Gamma: 0.0188   Vega: 37.5240
Parity check: C - P = 4.877058, S - K*exp(-rT) = 4.877058
Monte Carlo call: 10.4645 (+/- 0.0329), Black-Scholes: 10.4506
```

## Limitations
- European options only, no dividends
- Constant volatility and interest rate (no smile or skew)
