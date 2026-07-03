# ODE Validation: Monolithic vs Operator-Split

Compares the "ground truth" monolithic ODE solution for the dissolution model's
pressure diffusion against an approximate operator-split solver that mimics the
DGGML geocell decomposition.

## Model

11-node linear chain: `FluidSource -- Fluid -- ... -- Fluid -- FluidSink`

$$\frac{dP_i}{dt} = \alpha \sum_{j \in \text{neighbors}(i)} (P_j - P_i), \quad \alpha = 11.11$$

- Boundary: P₀ = 1.0 (fixed source), P₁₀ = 0.0 (fixed sink)
- Initial: P₀ = 1.0, all others = 0.0
- Steady state: P_i = 1 - i/10

## Usage

```bash
pip install numpy scipy matplotlib
python validate_ode.py
```

Produces `pressure_comparison.png` and prints max error at each time snapshot.
