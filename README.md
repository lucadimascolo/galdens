# galdens
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13928333.svg)](https://doi.org/10.5281/zenodo.13928333)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-31015/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![stability-alpha](https://img.shields.io/badge/stability-alpha-f4d03f.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#alpha)

A simple tool for generating galaxy density maps from a source list.

## Installation
> [!WARNING]
> This package was developed using Python 3.10. Newer versions of Python may work, but they have not been tested.


The package is not available on PyPI, so you will need to download it and install it locally. You can do this by cloning the repository and running the following command: 

```bash
git clone https://github.com/lucadimascolo/galdens.git
cd galdens
python -m pip install -e .
```

## Example
I am planning to work on a more extended documentation, but, for now, here is a basic example of how to use the `galdens` tool.

```python
import galdens
import numpy as np

from astropy import units as u

xpos = np.random.uniform(0,10) # RA
ypos = np.random.uniform(0,10) # Dec

# Use a histogram-based convolution approach
hdu_his = galdens.getmap(xpos,ypos,method='histogram'sigma=1*u.arcmin,cdelt=cdelt)

# Use Kernel Density Estimation with automatic kernel selection
hdu_kde = galdens.getmap(xpos,ypos,method='kde',cdelt=cdelt)

```
