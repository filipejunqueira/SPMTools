# SPMTools

SPMTools is a collection of scripts and utilities for processing and analyzing data from Scanning Probe Microscopy (SPM) experiments. The repository contains standalone programs for plotting spectroscopy curves, calculating forces from frequency shift measurements and batch processing SPM data files exported from various software packages.

Most of the scripts use data produced by ScientaOmicron's `Matrix` file format (`.mtrx`). The main helper functions for loading and analyzing these files are implemented in `spmUtils.py`. Example analysis results are stored in the repository as CSV files and images for reference.

## Repository layout

```
archive/               Miscellaneous scripts and older helper functions
curves/                Saved images of single spectroscopy curves
force_backup/          Backup of force calculations and plots
force_graphs/          Output plots of force calculations
several_curve_graphs/  Plots containing multiple curves
single_curve_graphs/   Plots of individual curves
*.csv                  Example data exported from the scripts
*.py                   Analysis and plotting scripts
```

A few key scripts are:

- **`spmUtils.py`** – Core library containing helper routines to load spectra, perform Savitzky–Golay filtering, fit curves, and calculate tip–sample forces using the Sader–Jarvis deconvolution algorithm.
- **`plot_single_curve.py`** – Example script showing how to average multiple spectra and plot a single curve with optional curve fitting and background subtraction.
- **`plot_several_curves.py`** – Plots multiple spectra in one figure using a color gradient for easy comparison.
- **`force_from_cvs.py`** – Computes short‑range forces from a list of ON/OFF spectroscopy pairs specified in CSV files.
- **`plot_trace_gwydion.py`** – Utility to plot line profiles exported from the Gwyddion SPM analysis package.

Example CSV lists such as `spec_list.csv` or `spec_on_off_list.csv` specify which files should be processed together.

## Requirements

The scripts rely on Python 3 with a number of scientific packages installed. At a minimum you will need:

- numpy
- pandas
- matplotlib
- seaborn
- scipy
- lmfit
- colour
- tkinter (for file dialogs)
- `access2thematrix` (third‑party module to read Omicron Matrix files)

Create a virtual environment and install dependencies with `pip`:

```bash
python3 -m venv venv
source venv/bin/activate
pip install numpy pandas matplotlib seaborn scipy lmfit colour
```

`access2thematrix` is not available on PyPI; obtain it from the original vendor or your internal distribution and ensure it is importable.

## Usage

The analysis workflow typically consists of:

1. **Preparing lists of spectra** – Edit `spec_list.csv` or the ON/OFF pair files to point to the spectra you want to process.
2. **Running the analysis scripts** – Execute one of the provided Python programs. For example:

   ```bash
   python plot_single_curve.py
   ```

   or to compute forces from ON/OFF pairs:

   ```bash
   python force_from_cvs.py
   ```

3. **Inspecting the results** – Output images will be saved in one of the folders mentioned above (for example `force_graphs` or `single_curve_graphs`). Many scripts also export processed data as CSV files for further inspection.

The scripts are intended to be edited to match your own directory layout and experiment parameters. Paths to data files as well as parameters such as Savitzky–Golay filter window, fitting algorithm and initial guesses are defined at the top of each script.

## Example: plotting a single curve

`plot_single_curve.py` demonstrates the typical processing steps:

1. Generate a list of file identifiers using `create_file_number_list` from `spmUtils.py`.
2. Average the curves with `average_curves`.
3. Optionally filter the data using `savgol_filter`.
4. Fit Lorentzian or Voigt peaks with `quick_fit`.
5. Plot the result with `plot_single_curve` and save the figure.

Parameters such as `curve_type`, `file_id`, number of curves, and the color scheme can be adjusted in the script.

## Notes

The repository contains many example CSV files produced during development. They serve as references for the expected output of the scripts but are not required for running new analyses.

Feel free to adapt the scripts to your own SPM experiments. Contributions and improvements are welcome!

