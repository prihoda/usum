# USUM: Plotting sequence similarity using USEARCH & UMAP

USUM uses [USEARCH](https://drive5.com/usearch/) and [UMAP](https://github.com/lmcinnes/umap) (or [t-SNE](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)) to plot DNA 🧬 and protein 🧶 sequence similarity embeddings.

[![PyPI - Downloads](https://img.shields.io/pypi/dm/usum.svg?color=green&label=PyPI%20downloads)](https://pypi.python.org/pypi/usum/)
[![PyPI license](https://img.shields.io/pypi/l/usum.svg)](https://pypi.python.org/pypi/usum/)
[![PyPI version](https://badge.fury.io/py/usum.svg)](https://pypi.python.org/pypi/usum/)
[![CI](https://api.travis-ci.org/prihoda/usum.svg?branch=master)](https://travis-ci.org/prihoda/usum)

## Installation

1. Install `USEARCH` dependency manually: https://drive5.com/usearch/download.html 
<br>(consider supporting the author by buying the 64bit license)

2. Install `usum` using PIP:

```bash
pip install usum
```

## Usage

Use `usum` to plot input protein or DNA sequences in FASTA format.

Show all available options using `usum --help`

### Minimal example


```bash
usum example.fa --maxdist 0.2 --termdist 0.3 --output example
```

### Multiple input files with labels

```bash
usum first.fa second.fa --labels First Second --maxdist 0.2 --termdist 0.3 --output example
```

This will produce a PNG plot:

![UMAP static example](docs/example1.png?raw=true "UMAP static example")

An interactive [Bokeh](https://bokeh.org) HTML plot is also created:

![UMAP Bokeh example](docs/example2.png?raw=true "UMAP Bokeh example")

### Using t-SNE instead of UMAP

You can also produce a t-SNE plot using the `--tsne` flag.

```bash
usum first.fa second.fa --labels First Second --maxdist 0.2 --termdist 0.3 --tsne --output example
```

This will produce a PNG plot:

![UMAP static example](docs/example1.png?raw=true "UMAP static example")
### Plotting random subset

You can use `--limit` to extract and plot a random subset of the input sequences.

```bash
# Plot 10k sequences from each input file
usum first.fa second.fa --labels First Second --limit 10000 --maxdist 0.2 --termdist 0.3 --output example
```

You can control randomness and reproducibility using the `--seed` option.

### Plotting options

See `usum --help` for all plotting options.

See [UMAP API Guide](https://umap-learn.readthedocs.io/en/latest/api.html) for more info about the UMAP options.

- Use `--limit` to plot a random subset of records
- Use `--width` and `--height` to control plot size in pixels
- Use `--resume` to reuse previous distance matrix from the output folder
- Use `--tsne` to produce a t-SNE embedding instead of UMAP (you can use this with `--resume`)
- Use `--umap-spread` to control how close together the embedded points are in the UMAP embedding
- Use `--umap-min-dist` to control minimum distance between points in UMAP embedding
- Use `--neighbors` to control number of neighbors in UMAP graph


### Reusing previous results

When changing just the plot options, you can use `--resume` to reuse previous results from the output folder.

**Warning** This will reuse the previous distance matrix, so changes to limits or USEARCH args won't take effect.

```bash
# Reuse result from umap output directory
usum --resume --output example --width 600 --height 600 --theme fire
```

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
