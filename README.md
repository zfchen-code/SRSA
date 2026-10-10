# SRSA supplementary code

MATLAB demonstrations accompanying the SRSA manuscript.

## Requirements

The code was developed using MATLAB R2024b Update 3 with Signal Processing Toolbox and Parallel Computing Toolbox. No GPU is required. `DoSRSA` uses `parfor`; parallel-pool settings follow the local MATLAB configuration. The reported computations used eight process workers.

## First run

Set the MATLAB **Current Folder** to the extracted `Supplementary_Code` directory, then run:

```matlab
addpath(fullfile(pwd,'src'));
cd('examples');
Demo2_FrequencyResolution
```

This adds `DoSRSA` to the MATLAB search path and runs Demo2. Run the other demos by entering their script names from the same `examples` directory. Keep `m1_2011.mat` in this directory for Demo4. Repeat the path setup in a new MATLAB session.

## Demonstrations

| Script | Purpose |
|---|---|
| `Demo1_KRobustness` | Sensitivity to the specified input model order. |
| `Demo2_FrequencyResolution` | Five- and nine-line synthetic examples corresponding to manuscript Figure 2(e, g). |
| `Demo3_NoiseResolution` | Nine groups of noisy two-tone signals, with 100 Monte Carlo realizations per group. |
| `Demo4_M1RealDataApplication` | Application to the supplied M1 superconducting-gravimeter record. |

Demo3 performs the full 900-realization experiment and takes substantially longer than the other demos. Runtime depends on hardware and parallel configuration; representative runtimes are given in the script headers.

## Files and use

- `src/`: the `DoSRSA` function, with complete SVD (`'svd'`, default) and randomized SVD (`'rsvd'`) modes.
- `examples/`: demonstration scripts and `m1_2011.mat`, the input data for Demo4.
- `figures/`: saved MATLAB figures.

Demo1, Demo2, and Demo4 use complete SVD; Demo3 uses randomized SVD for its SRSA calculations. The scripts generate workspace variables and display figures; they do not automatically overwrite the supplied FIG files.

The analysis parameters are set for the specified input series. If the sampling interval or input series is changed, reassess the input model order and spectral stacking interval, and update the corresponding row indices. See the script headers for details.
