# USUM: Plotting sequence similarity using USEARCH & UMAP

USUM uses [USEARCH](https://drive5.com/usearch/) and [UMAP](https://github.com/lmcinnes/umap) to plot DNA 🧬and protein 🧶 sequence similarity embeddings.

[![PyPI - Downloads](https://img.shields.io/pypi/dm/usum.svg?color=green&label=PyPI%20downloads)](https://pypi.python.org/pypi/usum/)
[![PyPI license](https://img.shields.io/pypi/l/usum.svg)](https://pypi.python.org/pypi/usum/)
[![PyPI version](https://badge.fury.io/py/usum.svg)](https://pypi.python.org/pypi/usum/)

## Installation

Install `UCLUST` manually: https://drive5.com/usearch/download.html (consider supporting the author by buying the 64bit license)

Install `usum` using PIP:

```bash
pip install usum
```

## Usage

### Minimal example

Plot input protein or DNA sequences in FASTA format:

```bash
usum sequences.fa --maxdist 0.2 --termdist 0.3 --output umap
```

### Multiple input files with labels

```bash
usum first.fa second.fa --labels First Second --maxdist 0.2 --termdist 0.3 --output umap
```

This will produce a PNG plot:

![UMAP static example](docs/example1.png?raw=true "UMAP static example")

An interactive [Bokeh](https://bokeh.org) HTML plot is also created:

![UMAP Bokeh example](docs/example2.png?raw=true "UMAP Bokeh example")

### Programmatic use

```python
from usum import usum

# Show help
help(usum)

# Run USUM
usum(inputs=['input.fa'], output='usum', maxdist=0.2, termdist=0.3)
```

## How it works

- A sparse distance matrix is calculated using USEARCH [calc_distmx](https://drive5.com/usearch/manual/cmd_calc_distmx.html) command. 
- The distances are based on % identity, so the method is agnostic to sequence type (DNA or protein)
- The distance matrix is embedded as a `precomputed` metric using [UMAP](https://github.com/lmcinnes/umap) 
- The embedding is plotted using [umap.plot](https://umap-learn.readthedocs.io/en/latest/plotting.html).
