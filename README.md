# respmethods - Robust circular cluster-based statistics for respiration-brain coupling
Robust pipeline for respiration phase-locked analyses of neural or behavioral data across the frequency spectrum, including a circular cluster permutation procedure reliably controlling error rates for multiple comparisons while accounting for the wrap-around continuity of circular data.

## Reference:
Berther, T., Balestrieri, E., Saltafossi, M., Paulsen, L. B., Andersen, L. M., Kluger, D. S. (2025). Robust circular cluster-based statistics for respiration-brain coupling. _PsyArXiv_, 2025-11.

Please cite this paper when the provided code is used.

## Contents:
- [MATLAB toolbox](matlab-toolbox---respLABmethods)
- [Python package](python-package---respymethods)

# MATLAB toolbox - respLABmethods
### Install
No special installation required, just add the 'matlab' folder in this repository and all its subfolders to your path.

### Functions
- ``CircClust`` cluster circular data points using complex representations in the polar plane
- ``CircPerm`` permutation test clusters in circular data
- ``two_point_interp``extract phase of an oscillatory time series using the two-point interpolation method
- ``four_point_interp``extract phase of an oscillatory time series using the four-point interpolation method
- ``protohase_interp`` extract phase of an oscillatory time series using the protophase interpolation method, adapted from the [DAMOCO toolbox](http://www.stat.physik.uni-potsdam.de/~mros/damoco.html) (C) Björn Kralemann, Michael Rosenblum, Arkady Pikovsky, University of Potsdam

#### Required, included functions
- ``generate_surrogate_iaaft``generates a phase-independent surrogate of a time series using the Iterative Amplitude-Adjusted Fourier Transform, (C) Alessio Perinelli & Leonardo Ricci, NSE Laboratory, Department of Physics, University of Trento ([github](https://github.com/LeonardoRicci/iaaft))
- ``PLV``computes the phase-locking value between two signals, (C) Edden Gerber, Edmond and Lily Safra Center for Brain Sciences, Hebrew University of Jerusalem, Israel ([github](https://github.com/edden-gerber/time_series_analysis_and_statistics))

### Scripts
- ``Tutorial_DataPrep.m`` demo pipeline for respiratory phase extraction, surrogate generation using IAAFT, and binning of behavioral data into respiratory phase bins
- ``Tutorial_DataClust.m`` demo script for circular permutation test on example empirical and simulated data, generates linear and polar plots of the results

### Example plots for results of circular cluster permutation test on randomly simulated data:

![results plot for simulated data, linear and polar representations](simcluster_results_matlab.png)

# Python package - respymethods
### Disclaimer
This package is still under active development and only available in a beta version, thus also does not yet have a Python Package Index. The code is not documented fully yet and partly experimental, download at your own risk.

### Install
#### 1. Requirements
The build of the package requires a C compiler, which may not necessarily be natively included in your OS. We recommend:
- Linux: `gcc`, install as recommended in the [gcc documentation](https://gcc.gnu.org/install/)
- MacOS: `clang`,  shipped with the [XCode Developer Tools](https://developer.apple.com/documentation/safari-developer-tools/installing-xcode-and-simulators)
- Windows:  `Visual C++`, download [here](https://visualstudio.microsoft.com/)


Additionally, the installation depends on the package manager uv - a fast and efficient python package and project manager. Please refer to the [documentation](https://docs.astral.sh/uv/) for more details.
On **Unix** systems (Linux + MacOS), the installer can be downloaded with `curl` and executed with `sh`:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

On **Windows**, use `irm` to download the executable script for `uv` and execute it with `iex` in a `PowerShell` prompt:

```bash
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

#### 2. Package installation
Respymethods will receive a Python Package Index in the future for easy installation using `pip`. The beta-dev version can be installed by manually compiling the included C modules and then installing the package. In line with best practices on python project management, we **highly** recommend creating your local installation of respymethods in a `uv` virtual environment.

Package installation works the same both on Unix and Windows platforms once `uv` has been installed.
From the parent repository folder (or the location where you usually stash your venvs), create a `uv` virtual environment in the `Terminal` (Unix) /`PowerShell` (Windows):

```bash
uv venv respyvenv
```

on **Unix**, activate it with:
```bash
source respyvenv/bin/activate
```
and on **Windows** with:
```powershell
respyvenv\Scripts\activate
```

Then, install the required packages in the virtual environment:

```bash
uv pip install setuptools numpy scipy matplotlib pyfftw
```

Manually compile the C modules by navigating into python/src/respymethods and running:

```bash
cd $PATH/respmethods/python/src/respymethods/
python setup.py build_ext --inplace
```

Upon successful compilation, `src/respymethods` then contains a `build` folder and `.so` files for the C modules.
Navigate back to the package root folder at `respmethods/python` and install using:

```bash
cd ../..
uv pip install -e .
```


### Functions
Detailed documentation of all included functions to be added soon.

### Scripts
- ``Tutorial_DataPrepMultiprocessing.py`` demo pipeline for parallelized, fast respiratory phase extraction and surrogate generation using IAAFT
- ``Tutorial_DataClust.py`` demo script for circular permutation test on example empirical and simulated data, generates linear and polar plots of the results


## Authors:
Teresa Berther, Elio Balestrieri, Daniel S. Kluger & Martina Saltafossi, Institute for Biomagnetism and Biosignal Analysis, University of Münster, Germany

Contact: teresa.berther@uni-muenster.de

## Disclaimer:
The code in this repository was implemented with care and tested on both empirical and simulated data. Nevertheless, it may contain errors or bugs, which may affect the outcome of your analysis. We do not take responsibility for any harm coming from using the provided code, neither caused by errors in the code nor by its improper application. Please email us any bugs you find.
