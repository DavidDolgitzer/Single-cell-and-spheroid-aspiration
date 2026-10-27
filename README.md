<h1 align="center">Contributions of single-cell mechanics and cell-cell adhesion to multicellular spheroid mechanics</h1>

<div align="center">

[David Dolgitzer](https://scholar.google.com/citations?user=SlyYieYAAAAJ&hl=en&oi=ao)<sup>1</sup>,
[Eleana Parajón](https://scholar.google.com/citations?user=y7IfxzUAAAAJ&hl=en&oi=ao)<sup>2</sup>,
[Douglas N. Robinson](https://scholar.google.com/citations?user=vuHsbR8AAAAJ&hl=en&oi=ao)<sup>2</sup>,
[Pablo A. Iglesias](https://scholar.google.com/citations?user=iRdtkAgAAAAJ&hl=en)<sup>3</sup>

<sup>1</sup> Department of Physics and Astronomy, Krieger School of Arts and Sciences, Johns Hopkins University, Baltimore, MD, 21218, USA  
<sup>2</sup> Department of Cell Biology, School of Medicine, Johns Hopkins University, Baltimore, MD, 21205, USA  
<sup>3</sup> Department of Electrical and Computer Engineering, Whiting School of Engineering, Johns Hopkins University, Baltimore, MD, 21218, USA

</div>

This repository contains simulation models, simulation outputs, analysis code, and plotting notebooks for the study. It supports three levels of reproduction:

1. **Recreate the manuscript figures from the included `.mat` files.** Morpheus and MATLAB are not required for this route.
2. **Repeat the analysis from the included CSV simulation outputs.** This route requires MATLAB.
3. **Run new simulations from the included Morpheus XML models.** This route requires Morpheus, followed by MATLAB and Python for analysis and plotting.

## Repository contents

- `Single-cell aspiration/` — single-cell simulation outputs, MATLAB analysis, and plotting notebook.
- `Spheroid aspiration/` — spheroid simulation outputs, MATLAB analysis, and plotting notebooks.

Keep the folder structure unchanged: the MATLAB scripts and Python notebooks use relative paths.

## Simulation conditions

The spheroid analysis covers 25 combinations of single-cell stiffness and cell-cell adhesion:

- `λₛ`: 0.8, 0.9, 1.0, 1.1, and 1.2
- `J`: -40, -45, -50, -55, and -60

Here, `λₛ` represents [FILL IN: meaning and units, if applicable], and `J` represents [FILL IN: meaning and units, if applicable].

Folder names use this format:

    lambda_S_[stiffness]_J_[adhesion]

For example, `lambda_S_1_J_50` contains simulations with `λₛ = 1` and `J = 50`.

## Recreate figures from the processed data

The processed `.mat` files are included so that figures can be recreated without Morpheus or MATLAB.

### Requirements

Tested environment:
- Python (Jupyter Notebook): 3.11.5 
- NumPy: 1.26.4
- SciPy: 1.11.1
- Matplotlib: 3.7.2
- MATLAB: R2023a Update 5 (9.14.0.2337262)

### Figure 2: single-cell aspiration

Open and run:

    Single-cell aspiration/Python figure/Figure2_single_cell_elasticity.ipynb

The notebook reads processed `.mat` files from these single-cell stiffness conditions:

    lambda_S_0.8
    lambda_S_1
    lambda_S_1.2

### Figures 4 and 5: spheroid aspiration

Open and run:

    Spheroid aspiration/Python figure/Figure4_single_cell_contribution.ipynb
    Spheroid aspiration/Python figure/Figure5_adhesion_contribution.ipynb

These figures use five parameter combinations:

    lambda_S_1_J_50
    lambda_S_1_J_40
    lambda_S_1_J_60
    lambda_S_0.8_J_50
    lambda_S_1.2_J_50

These combinations contain the additional aspiration-force conditions required for Figures 4 and 5.

### Figure 6: parameter-space heatmaps

Open and run:

    Spheroid aspiration/Python figure/Figure6_tau_heatmaps.ipynb

This notebook uses all 25 combinations of `λₛ` and `J`, with aspiration-force conditions, to calculate and display retardation-time (`τ`) heatmaps.

Run each notebook from beginning to end. The notebooks use repository-relative paths, so keep the directory structure unchanged.

## Repeat the MATLAB analysis from CSV files

### Requirements

- MATLAB
- Optimization Toolbox (`lsqcurvefit`)
- Statistics and Machine Learning Toolbox (`nlparci` and `tinv`)

### Single-cell aspiration analysis

Each `force_0.XX` directory contains outputs for 10 cells simulated simultaneously. The MATLAB analysis aligns the 10 aspiration-length curves, calculates their mean and standard error, fits a Kelvin–Voigt model, and estimates mechanical parameters and confidence intervals.

To repeat the analysis for one stiffness condition:

1. Open its `sim_analysis_MATLAB` directory in MATLAB.
2. Run `process_all_blocks.m`.
3. Confirm that the corresponding `block_XX_results.mat` files are created or updated.

The `force_XX.m` scripts connect each aspiration-force directory to `process_sim_block.m`. The resulting `.mat` files are used by the Python plotting notebook.

Script must be run from its own directory. Changes in analysis code will change the existing `block_XX_results.mat`, therefore changing the figure plots in python.]

### Spheroid aspiration analysis

Each `sim#` directory represents one independent Morpheus simulation. Spheroid simulations were run independently because they were substantially more computationally expensive than the single-cell simulations. At least 10 independent simulations contribute to each force condition.

The spheroid MATLAB analysis aligns the independent simulations, applies the boundary-handling procedure where required, calculates ensemble statistics, fits the spheroid viscoelastic model, and saves fitted parameters and plotting data in `.mat` files.

To repeat the analysis for one parameter combination:

1. Open that combination’s `sim_analysis_MATLAB` directory in MATLAB.
2. Run `process_all_blocks.m`.
3. Confirm that the corresponding `block_[force]_results.mat` files are created or updated.

Script must be run from its own directory. Changes in analysis code will change the existing `block_XX_results.mat`, therefore changing the figure plots in python.]

## Run simulations in Morpheus

### Requirements

- Morpheus Version 2.3.9
- 4 threads with random seed

Two Morpheus XML models are provided:

1. **Single-cell aspiration model**
   - Default stiffness: `λₛ = 1`
   - Change the `SurfaceConstraint` parameter to reproduce other single-cell stiffness conditions, including 0.8 and 1.2.

2. **Spheroid aspiration model**
   - Default single-cell stiffness: `λₛ = 1`
   - Default cell-cell adhesion: `J = -50`
   - Change the corresponding stiffness and cell-cell contact-adhesion parameters to explore other combinations.

To run a model:

1. Install Morpheus.
2. Open the relevant XML file.
3. Set the desired stiffness, adhesion, aspiration force, and random seed.
4. Run the simulation.
5. Save the generated CSV output in the corresponding directory structure.
6. Run the MATLAB analysis, then the appropriate Python notebook.

Parameter names in the XML files:

- Aspiration force: `DirectedMotion strength`
- `λₛ`:   `SurfaceConstraint strength`
- `J`: `cell-cell Contact`

## Data and output conventions

- CSV files are raw Morpheus simulation outputs.
- Each single-cell force condition contains 10 cell-level CSV files from one simultaneous simulation.
- Each spheroid `sim#` folder represents one independent simulation run and contains one CSV output.
- `.mat` files are processed analysis outputs generated by MATLAB and included so that the Python plotting workflow can be used without MATLAB.
- Python notebooks generate manuscript plots from the processed `.mat` files.
- Force labels in filenames: for example `force_06` = 0.6 and `force_1` = 1.0.

## Reproducibility notes

- Keep the repository directory structure unchanged because the MATLAB scripts and Python notebooks use relative paths.
- Use 4 threads and keep seed random (the latter is a default in Morpheus).
- Results may vary when simulations are rerun with different random seeds.
- The included `.mat` files were generated from the CSV outputs included in this repository.
- Approximate runtime and hardware: single-cell simulations: 15 seconds; 2:30-5 hours, depending on aspiration waiting time. These times represent simulations run on a 2021 MacBook Pro with an Apple M1 processor and 16 GB RAM.

### Manually assigned parameters

Several parameters were manually assigned during the `spheroid simulation` analysis in the MATLAB scripts
- `clipping` (`force_#.m`): clips aspiration length and establishes t=0 from there. This was done manually to ignore initial lags in aspiration. Generally, lower aspiration force, stiffer cells, and stronger adhesion require more clipping.
- `analysis_time_limit` (`force_#.m`): tracks the slowest simulation to reach the end of the pipette, and extrapolates all other simulation to reach that time. Generally, weaker adhesion, softer cells, and higher aspiration force requires extrapolation.
- `baseline_force` (`process_sim_block`): The amount of aspiration required for initial aspiration, sometimes called "critical pressure". To find, plot `mean_L_end` vs. ΔF (as in Figure 4d, only with `mean_L_end` instead of `mean_L`) and fit to find the y-axis intersection. Preferably use only lower aspiration force.

## License

Code is licensed under the MIT License. See LICENSE.
Morpheus is a separate software project; see its official website and license for its terms.

## Authors and contact

For more information discussion please contact David Dolgitzer at ddolgit1@jhu.edu.
