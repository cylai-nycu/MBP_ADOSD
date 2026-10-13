# DS-CSS MBP4 + ADOSD decoder

This is a Windows command-line decoder built from the upstream C implementation in `ADOSD/ldpcQaun`. It keeps the original MBP4 message passing, ADOSD post-decoder, reliability sorting, and RSR search. The old simulation modes are not part of this dedicated executable.

Supported inputs:

- `18-original`: `[[18,6,3]]` with the original six X and six Z generators.
- `18-ds3`: `[[18,6,3]]` with the DS-distance-3 generator set.
- `36`: `[[36,12,5]]` with sixteen X and sixteen Z generators.
- `36-4-6`: `[[36,4,6]]` with 23 X and 23 Z measurements, including the supplied sparse redundant checks (DS-distance 5).

## Build on Windows

Install MinGW-w64 GCC and ensure `gcc.exe` is on `PATH`. From the repository root run:

```powershell
.\DS_CSS\ldpcQaun\build.bat
```

Or, if already in `DS_CSS\ldpcQaun`:

```powershell
.\build.bat
```

The build creates `DS_CSS\ldpcQaun\DS_CSS_decoder_parallel.exe`. Run it from that directory so it can load the bundled matrices.

## Run

```powershell
cd .\DS_CSS\ldpcQaun
.\DS_CSS_decoder_parallel.exe --code 18-ds3 --p-data 0.01 --p-meas 0.01 --iterations 50 --order 4 --shots 1000
.\DS_CSS_decoder_parallel.exe --code 36 --p-data 0.01 --p-meas-x 0.008 --p-meas-z 0.012 --iterations 50 --order -2 --shots 1000
.\DS_CSS_decoder_parallel.exe --code 36-4-6 --p-data 0.0001 --p-meas 0.0001 --iterations 100 --order 2 --target-failures 400 --shots 2147483647 --threads 8
```

Options:

- `--code`: `18-original`, `18-ds3`, `36`, or `36-4-6` (required).
- `--p-data`: total physical depolarizing probability. Physical X, Y, and Z each have probability `p_data/3`.
- `--p-meas`: default probability for both virtual sectors.
- `--p-meas-x`: optional override for X-syndrome virtual qubits, which allow only Z errors.
- `--p-meas-z`: optional override for Z-syndrome virtual qubits, which allow only X errors.
- `--iterations`: MBP iteration limit, 1–100.
- `--order`: `-2` selects adaptive APOSD with RSR reduction, `-1` disables ADOSD, and `0..12` select fixed OSD order. The default is `4` for `18-ds3` and `2` for the other codes. Fixed orders use the upstream `post_decodeOSDw` reliability-sorted path; `-2` uses `post_decodeOSDfull` with RSR.
- `--shots`: Monte Carlo shots.
- `--threads`: Monte Carlo worker threads (default 1). Each worker owns separate BP and ADOSD buffers and a separate random stream; results are combined after all workers finish. Set this to the number of logical CPU cores you want to use.
- `--target-failures`: stop once this many logical failures are observed, up to the `--shots` maximum. The output reports the actual number of shots and whether the target was reached.
- `--always-osd`: always invoke ADOSD after BP, including when BP's hard decision already matches the syndrome. Without this option, a syndrome-consistent BP result is returned directly. This option cannot be combined with `--order -1`.
- `--seed`: reproducible channel-sampling seed.
- `--self-test`: check syndrome generation, single faults, physical normalizer classification, and zero-probability virtual Pauli components. DS3 and `[[36,12,5]]` also run end-to-end decoding on every single physical and measurement fault.

Output includes `logical_error_rate`, `failures`, `shots`, and `average_iterations`. A failure is counted only when the data-qubit residual anticommutes with a physical normalizer row. Virtual measurement correction differences do not count as logical errors.

Every Monte Carlo run also appends one row to a text file named after the selected code (`18-original.txt`, `18-ds3.txt`, `36.txt`, or `36-4-6.txt`) in the current working directory. The file starts with `//p_data\tlogical_error_rate\terror_bar`; each subsequent row contains those three values. `error_bar` is the approximate 95% binomial interval, calculated as `2 * sqrt(LER * (1 - LER) / shots)` (two standard errors). Repeated runs append rows without replacing earlier results.

## Decoder and channel details

The matrices are bundled in `DS_CSS/codes/18-original`, `DS_CSS/codes/18-ds3`, `DS_CSS/codes/36`, and `DS_CSS/codes/36-4-6`. `DS_A4.txt` is the augmented parity-check matrix. Its variables are physical data qubits, X-syndrome virtual qubits, then Z-syndrome virtual qubits. The X-syndrome virtual sites have prior `{I:1-p_meas_x, Z:p_meas_x}`; the Z-syndrome sites have `{I:1-p_meas_z, X:p_meas_z}`. Forbidden Pauli probabilities stay exactly zero. Dummy capacity added for the fixed-size upstream arrays has identity probability 1 and no check edges.

The original BP/OSD core is compiled with storage for the largest supported input (`N=82`, `M=46`, `M_G=48`); BP and OSD operate on each code's actual runtime dimensions. After decoding, virtual corrections are recomputed from the observed syndrome and physical estimate, then checked for full syndrome consistency and legal Pauli support. The physical normalizer in `DS_G4s.txt` has zero virtual columns, so logical success is measured only on data qubits.

Adaptive `-2` APOSD can take substantially longer on difficult syndromes. Fixed order 4 is the default for `18-ds3` to ensure the supported DS single-fault suite is corrected at the tested priors. Run the bundled checks from the repository root with:

```powershell
powershell -ExecutionPolicy Bypass -File .\DS_CSS\run_tests.ps1
```

The decoder sources are distributed under the upstream repository's GPL-3.0 license; see `LICENSE`.
