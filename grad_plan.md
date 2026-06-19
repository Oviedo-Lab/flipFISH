# Plan: Analytical Gradient of `mQC_msle` for nlopt

## Computation Chain

`mQC_msle(params, data_ptr)` computes the following pipeline:

```
params (vector of length 2*N_bits + 2*corr_free)
  → pack_fr()             → FlipRates fr
  → est_bc_counts_true()  → bc_counts_true   [integer vector — treat as constant]
  → expected_bc_counts()  → ecc              [vector of predicted corrected counts]
  → compute_msle()        → scalar MSLE
```

Parameter layout (as packed by `pack_fr`):
| Index range                        | Parameter               |
|------------------------------------|-------------------------|
| `[0, N_bits)`                      | `rate10[i]`  (P(1→0))   |
| `[N_bits, 2*N_bits)`               | `rate01[i]`  (P(0→1))   |
| `[2*N_bits, 2*N_bits+corr_free)`   | `corr1[k]`   (k = i*(i-1)/2+j) |
| `[2*N_bits+corr_free, end)`        | `corr0[k]`              |

---

## Mathematical Derivation

### MSLE gradient (outer layer)

With uniform weights and observed counts `c_b` (fixed):

$$\frac{\partial \text{MSLE}}{\partial \theta_k}
  = \frac{1}{N_B} \sum_b
      \frac{2\bigl[\log(\hat{c}_b+1) - \log(c_b+1)\bigr]}{\hat{c}_b + 1}
      \cdot \frac{\partial \hat{c}_b}{\partial \theta_k}$$

### Expected corrected count gradient (middle layer)

`expected_bc_count` for barcode `b` sums over all codebook barcodes `j` and all raw reads `i` that are error-corrected to `b`:

$$\hat{c}_b = \sum_j n_j^{\text{true}} \sum_{i:\,\text{corr}(i)=b} TR(j,\; j \oplus i,\; \theta)$$

so:

$$\frac{\partial \hat{c}_b}{\partial \theta_k}
  = \sum_j n_j^{\text{true}} \sum_{i:\,\text{corr}(i)=b} \frac{\partial TR(j, j \oplus i, \theta)}{\partial \theta_k}$$

**Key decision:** Treat `bc_counts_true` as a constant w.r.t. differentiation. It is computed via rounding (a step function), so it is non-differentiable anyway, and its sensitivity to params is weak compared to the TR terms.

### TR gradient (inner layer — the core derivation)

`TR(bc, flips, fr)` computes:

$$TR = \prod_i f_i, \quad
  f_i = \begin{cases} p_i & \text{if } \text{flip}_i = 1 \\ 1 - p_i & \text{if } \text{flip}_i = 0 \end{cases}$$

where the **correlation-adjusted flip probability** for bit `i` is:

$$p_i = \exp\!\bigl(\underbrace{\log\text{rate\_flip}_i}_{\text{log\_rate10}[i] \text{ or } \text{log\_rate01}[i]} - \underbrace{\sum_{j < i,\, \text{cond}(i,j)} \log(1 - \text{corr}_{ij})}_{\text{log\_inv\_corr}[i]}\bigr)$$

where the condition `cond(i,j)` is:
- use `corr1[k]` when `bit_i = 1` AND `flip_j = 1`
- use `corr0[k]` when `bit_i = 0` AND `flip_j = 1`

Since `log TR = Σ_i log f_i`, the gradient decomposes bit-by-bit. Define the **sensitivity factor** for bit `i`:

$$\delta_i \;=\; \frac{\partial \log TR}{\partial \log p_i}
  = \begin{cases} 1 & \text{if } \text{flip}_i = 1 \\ -p_i/(1-p_i) & \text{if } \text{flip}_i = 0 \end{cases}$$

Then: $\partial TR / \partial \theta = TR \cdot \partial \log TR / \partial \theta$.

#### Gradient w.r.t. `rate10[i]` (only relevant when `bit_i = 1`):

$$\frac{\partial TR}{\partial \text{rate10}_i}
  = TR \cdot \delta_i \cdot \frac{1}{\text{rate10}_i}$$

#### Gradient w.r.t. `rate01[i]` (only relevant when `bit_i = 0`):

$$\frac{\partial TR}{\partial \text{rate01}_i}
  = TR \cdot \delta_i \cdot \frac{1}{\text{rate01}_i}$$

#### Gradient w.r.t. `corr1[k]` where `k = i*(i-1)/2 + j`:

Non-zero only when `bit_i = 1` AND `flip_j = 1` (i.e., when `corr1[k]` appears in `log_inv_corr[i]` for this `(bc, flips)` pair).

$$\frac{\partial TR}{\partial \text{corr1}_k}
  = TR \cdot \delta_i \cdot \frac{1}{1 - \text{corr1}_k}$$

#### Gradient w.r.t. `corr0[k]`:

Non-zero only when `bit_i = 0` AND `flip_j = 1`:

$$\frac{\partial TR}{\partial \text{corr0}_k}
  = TR \cdot \delta_i \cdot \frac{1}{1 - \text{corr0}_k}$$

#### Clamping region:

When `log_adj_flip[i] >= 0` (p_i would exceed 1), the existing code clamps to a constant. In the gradient, treat this as a **zero derivative** for all parameters affecting bit `i`'s log_adj_flip (i.e., set `δ_i = 0`).

---

## New Functions to Write

### 1. `TR_and_grad`

```cpp
// Computes TR(bc, flips, fr) and accumulates ∂TR/∂params into grad_* vectors.
// The grad_* vectors are ADDITIVE — initialize to zero before the first call.
double TR_and_grad(
    uint64_t bc,
    uint64_t flips,
    const FlipRates& fr,
    int N_bits,
    std::vector<double>& grad_rate10,   // length N_bits
    std::vector<double>& grad_rate01,   // length N_bits
    std::vector<double>& grad_corr1,    // length corr_free
    std::vector<double>& grad_corr0     // length corr_free
);
```

Implementation mirrors `TR` exactly for the forward pass, then adds an O(N_bits²) backward pass to fill the gradient vectors.

### 2. `expected_bc_count_and_grad`

```cpp
// Returns (count_read, count_corrected, count_true) AND fills grad_ecc (∂ecc[b]/∂params).
std::tuple<double, double, double> expected_bc_count_and_grad(
    uint64_t BOI,
    const FlipRates& fr,
    const std::vector<int>& bc_counts_true,
    const std::vector<uint64_t>& true_barcodes,
    const std::vector<uint64_t>& corrected_to_BOI,
    std::vector<double>& grad_ecc        // length n_params, ADDITIVE
);
```

Internally calls `TR_and_grad` instead of `TR` for each (bc_true, flips) pair.

### 3. `expected_bc_counts_and_grad`

Mirrors `expected_bc_counts` but also returns a matrix of gradients
`decc[b][k]` = ∂ecc[b]/∂θ_k. If parallelism (n_forks > 1) is needed, each child process can write back an extended buffer containing both the three count vectors and the full gradient matrix for its batch of barcodes.

### 4. `mQC_msle_grad` (nlopt-compatible)

```cpp
// [[Rcpp::export]]
double mQC_msle_grad(
    NumericVector params,
    SEXP data_ptr_sexp,
    bool compute_gradient,           // if false, behaves like mQC_msle
    NumericVector grad_out           // modified in-place when compute_gradient=true
);
```

Or alternatively, a free function with the C-style nlopt callback signature:
```cpp
double nlopt_objective(
    unsigned n,
    const double* x,
    double* grad,       // NULL if gradient not needed
    void* f_data
);
```

The body:
1. `pack_fr(params)` → `fr`
2. `est_bc_counts_true(fr, d)` → `bc_counts_true`  ← treated as constant
3. `expected_bc_counts_and_grad(fr, bc_counts_true, ...) ` → `(ecc, decc)`
4. Apply MSLE chain rule to combine `ecc`, `d->bc_counts`, and `decc` into `grad_out`
5. Return MSLE scalar

---

## Order of Implementation

1. **`TR_and_grad`** — self-contained, easy to unit-test against finite differences on small examples.
2. **`expected_bc_count_and_grad`** — wraps `TR_and_grad`, testable the same way.
3. **`expected_bc_counts_and_grad`** — decide upfront whether to include parallel (fork) support or serial-only first; serial is simpler and still useful for nlopt.
4. **`mQC_msle_grad`** — final assembly; export with Rcpp.

---

## Testing Strategy

For each new function, verify against finite-difference numerical gradients:

```r
fd_grad <- function(f, x, h = 1e-6) {
  sapply(seq_along(x), function(i) {
    xp <- x; xp[i] <- xp[i] + h
    xm <- x; xm[i] <- xm[i] - h
    (f(xp) - f(xm)) / (2 * h)
  })
}
```

Use a small synthetic codebook (N_bits = 4, a handful of barcodes) to keep the test fast.

---

## Open Questions / Caveats

- **`bc_counts_true` discontinuity:** The rounding in `est_bc_counts_true` makes MSLE technically non-differentiable at parameter values where integer rounding changes. In practice this is rare and the gradient from treating true counts as fixed should work well for nlopt's gradient-based methods (e.g., L-BFGS). Worth noting in comments.
- **Corr parameter range:** `corr1[k]` and `corr0[k]` must stay in `(-∞, 1)` for `log(1 - corr)` to be real. nlopt bounds should enforce this; the gradient formula has a singularity at `corr = 1`.
- **Parallel fork support in gradient:** Passing back full gradient matrices through pipes is straightforward but increases pipe buffer size by a factor of `n_params`. For large N_bits this may exceed default pipe buffer limits; may need to loop reads/writes (already done in the existing pipe code).