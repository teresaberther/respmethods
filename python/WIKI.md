# respymethods

## How to get started

The package is still in a pre-beta version, so it does not have yet a Python Package Index. 
The following guide for using respymethods has only been tested on Linux.

### Install uv

uv is an extremely fast Python package and project manager. Please refer to the [documentation](https://docs.astral.sh/uv/) for more details.

On macOS and Linux you can download the installer with `curl` and execute it with `sh`:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Create the environment

You can create a local installation of respymethods in a virtual environment calling `uv sync` from the parent repo folder: 

```bash
cd respmethods/python
uv sync
```
And activate it with:

```bash
source .venv/bin/activate
```

## How to use respymethods

Check out `Tutorial_CircClust.py` and `Tutorial_CircPermMultiprocessing.py`.

