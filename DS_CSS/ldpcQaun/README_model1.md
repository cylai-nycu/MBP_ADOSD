# DS-CSS MBP4 + ADOSD: weight-dependent syndrome noise (Model 1)

This is a separate executable and entry point. It leaves `DS_CSS.c`, the shared BP/ADOSD implementation, and the supplied matrices untouched.

## Build on Windows

Install MinGW-w64 GCC with OpenMP and make sure `gcc.exe` is on `PATH`. From this directory, run:

```powershell
.\build_model1.bat
```

Run the executable from this directory so it can find `..\codes\...`:

```powershell
.\DS_CSS_model1.exe --code 18-ds3 --p-data 0.001 --iterations 50 --order 2 --shots 100000
.\DS_CSS_model1.exe --code 36 --p-data 0.001 --iterations 50 --order 2 --shots 100000 --threads 8
```

Supported code names are `18-original`, `18-ds3`, `36`, and `36-4-6`.

## Model 1

Edit `#define MEAS_RATIO 1.0` near the top of `DS_CSS_model1.c` to choose the measurement-noise ratio. The executable sets `p_m = MEAS_RATIO * p_data` and derives a separate probability for each check from the physical columns of the loaded augmented matrix:

```text
q_i = (1 - (1 - 2*p_m)^w_i) / 2
```

Here `w_i` counts nonzero entries in check `i` among the physical data-qubit columns only; virtual columns are excluded. The code rejects values where `p_m` is outside `[0, 0.5]`. Each virtual measurement bit is sampled independently with its own `q_i`; the same probability is installed in BP's variable priors and edge reliabilities, then passed to ADOSD through its posterior probabilities. X-syndrome virtual qubits allow only Z errors and Z-syndrome virtual qubits allow only X errors. Forbidden Pauli probabilities remain exactly zero.

The physical channel remains depolarizing: `P(I)=1-p_data`, `P(X)=P(Y)=P(Z)=p_data/3`. BP's existing syndrome-consistent early return and ADOSD fallback remain in place. The existing physical logical-error checker is retained.

## Options and output

- `--code NAME`: required; `18-original`, `18-ds3`, `36`, or `36-4-6`.
- `--p-data P`: physical depolarizing probability; measurement probability is derived from `MEAS_RATIO`.
- `--iterations N`: MBP iteration limit.
- `--order N`: ADOSD order (`-2` adaptive APOSD/RSR, `-1` disabled, or fixed order `0..12`).
- `--shots N`: maximum Monte Carlo samples.
- `--target-failures N`: optional stopping target, limited by `--shots`.
- `--threads N`: OpenMP workers.
- `--seed N`: random seed.
- `--self-test`: matrix/single-fault consistency, logical-classification, and virtual-prior checks; supported codes also run single-fault decoding checks.

The executable prints logical error rate, failure count, shots, and average BP iterations. Each simulation appends `p_data`, logical error rate, and an approximate two-standard-error binomial bar to `<code>.txt` in the current directory.
