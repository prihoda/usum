#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import shutil
import sys
from scipy import sparse
import pandas as pd
import umap
import bokeh.plotting
from Bio import SeqIO
from collections import OrderedDict
import numpy as np
import umap.plot
import warnings
import itertools
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt
from umap.plot import _matplotlib_points, _themes, _select_font_color, _datashade_points

def main(argv=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", nargs='*', help="Input FASTA sequence file paths.")
    parser.add_argument("--labels", nargs='*', required=False, help="Input file labels.")
    parser.add_argument("--output", required=True, help="Output directory path.")
    parser.add_argument("-f", "--force", action="store_true", default=False, help="Force overwrite output.")
    parser.add_argument("--resume", action='store_true', default=False, help="Resume using existing results from output directory.")

    parser.add_argument("--limit", type=int, help="Use random number of records from each input file.")
    parser.add_argument("--seed", default=1, type=int, help="Random seed for input subsampling and UMAP.")
    
    parser.add_argument("--maxdist", type=float, help="USEARCH: Maximum distance which should be written (required if not using --resume).")
    parser.add_argument("--termdist", type=float, default=1.0, help="USEARCH: Identity threshold for terminating the calculation. This should be set higher than maxdist.")
    
    parser.add_argument("--umap-min-dist", type=float, default=0.1, help="UMAP: Effective minimum distance between embedded points, relative to spread.")
    parser.add_argument("--umap-spread", type=float, default=1.0, help="UMAP: The effective scale of embedded points.")
    parser.add_argument("--neighbors", type=int, default=15, help="UMAP: The size of local neighborhood.")

    parser.add_argument("--theme", default='fire', help="Plot color theme.")
    parser.add_argument("--width", type=int, default=800, help="Plot width in pixels.")
    parser.add_argument("--height", type=int, default=800, help="Plot height in pixels.")
    
    parser.add_argument("--tsne", action="store_true", default=False, help="Run t-SNE instead of UMAP.")

    options = parser.parse_args()
    
    warnings.filterwarnings("ignore", message="using precomputed metric")

    if not options.resume and not options.inputs:
        parser.error("Input file paths are required when not using --resume")

    try:
        usum(
            inputs=options.inputs, 
            output=options.output, 
            maxdist=options.maxdist, 
            termdist=options.termdist,
            labels=options.labels, 
            force=options.force,
            resume=options.resume,
            limit=options.limit,
            random_state=options.seed,
            umap_min_dist=options.umap_min_dist,
            umap_spread=options.umap_spread,
            neighbors=options.neighbors,
            method='tsne' if options.tsne else 'umap',
            theme=options.theme,
            width=options.width,
            height=options.height
        )
    except UsumError as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)

class UsumError(Exception):
    pass
    
def usum(
        inputs, output, maxdist=None, termdist=1.0, 
        labels=None, force=False, resume=False, limit=None, random_state=1,
        umap_min_dist=0.1, umap_spread=1.0, method='umap', neighbors=15, theme='fire', width=800, height=800
    ):
    """
    Compute sequence similarity and plot UMAP embedding.
    :param inputs: list of FASTA input paths.
    :param output: output directory path.
    :param maxdist: USEARCH: Maximum distance which should be written.
    :param termdist: USEARCH: Identity threshold for terminating the calculation. This should be set higher than maxdist.
    :param labels: Input file labels. If not provided, file names without extension will be used.
    :param force: Force overwrite output.
    :param resume: Resume using existing results from output directory.
    :param limit: Use random number of records from each input file.
    :param random_state: Random seed for input subsampling and UMAP.
    :param umap_min_dist: UMAP: Effective minimum distance between embedded points
    :param umap_spread: UMAP: The effective scale of embedded points
    :param neighbors: UMAP: The size of local neighborhood.
    :param method: Embedding method (umap, tsne).
    :param theme: Plot color theme.
    :param width: Plot width in pixels.
    :param height: Plot height in pixels.
    :return: tuple with UMAP reducer and sequence DataFrame
    """
    if not force and not resume and (os.path.exists(output) and (os.listdir(output) or not os.path.isdir(output))):
        raise UsumError(f'Output path exists and is not empty: {output}. \nRun with -f --force to overwrite or --resume to resume')

    if labels:
        if len(labels) != len(inputs):
            raise ValueError('Labels have to be provided for each input file (--labels a b c)')
    else:
        labels = [os.path.splitext(os.path.basename(path))[0] for path in inputs]

    fasta_path = os.path.join(output, 'input.fa')
    distance_path = os.path.join(output, 'distance.txt')
    index_path = os.path.join(output, 'sequences.tsv')
    png_path = os.path.join(output, method+'.png')
    html_path = os.path.join(output, method+'.html')
        
    if resume and os.path.exists(index_path) and os.path.exists(distance_path):
        print(f'> Resuming using previous results...')
        if limit is not None:
            print('WARNING: --limit specified but is ignored')
        if maxdist is not None:
            print('WARNING: --maxdist specified but is ignored')
        index = pd.read_csv(index_path, sep='\t')
    else:
        print(f'> Writing FASTA records from {len(inputs)} paths...')
        if resume:
            print('WARNING: No previous results found, computing from scratch!')

        if not inputs:
            raise UsumError('Input file paths are required')
            
        if maxdist is None:
            raise UsumError('Argument --maxdist is required')
                    
        if not os.path.exists(output):
            os.mkdir(output)

        if os.path.exists(index_path):
            # Remove sequence TSV file to avoid reusing incomplete result
            os.remove(index_path)
            
        index = save_input_fasta(inputs, labels, fasta_path, limit=limit, random=(limit is not None), random_state=random_state)
        print(f'Saved {len(index):,} records to: {fasta_path}')

        print(f'\n> Creating sparse {len(index):,} x {len(index):,} distance matrix with {maxdist} max distance...')
        if len(index) > 10000 and not limit:
            print('NOTE: This might take some time. Consider using --limit to compare just a random subset.')
        run_usearch(fasta_path, distance_path, maxdist=maxdist, termdist=termdist)
    
    dist_matrix = load_sparse_dist_matrix(distance_path)
    # with {dist_matrix.getnnx():,} non-null values
    print(f'Loaded {len(index):,} x {len(index):,} distance matrix ({int(dist_matrix.nbytes / 1024 / 1024)} MB)')
    
    if method == 'tsne':
        print(f'\n> Creating t-SNE embedding...')
        reducer, embedding = fit_tsne(dist_matrix, random_state=random_state)
                
        index['tsne1'] = embedding[:,0]
        index['tsne2'] = embedding[:,1]
    elif method == 'umap':
        print(f'\n> Creating UMAP embedding with {neighbors} neighbors...')
        reducer, embedding = fit_umap(dist_matrix, neighbors=neighbors, random_state=random_state, min_dist=umap_min_dist, spread=umap_spread)
                
        index['umap1'] = embedding[:,0]
        index['umap2'] = embedding[:,1]
    else:
        raise ValueError(f'Unknown embedding method: {method}')

    index.to_csv(index_path, sep='\t', index=False)
    print(f'Saved sequences TSV to: {index_path}')

    print('\n> Drawing PNG...')
    ax = plot_points(reducer, labels=index['label'], theme=theme, width=width, height=height);
    ax.figure.savefig(png_path, bbox_inches='tight')
    print(f'Saved PNG to: {png_path}')
    
    print('\n> Drawing interactive plot...')
    p = umap.plot.interactive(reducer, labels=index['label'], theme=theme, width=width, height=height, hover_data=index);
    bokeh.plotting.output_file(html_path)
    bokeh.plotting.save(p)
    print(f'Saved plot HTML to: {html_path}')
    
    print(f'\nDone. Saved to: {output}')
    return reducer, index

def iterate_fasta_index(records, idx):
    for i, record in enumerate(records):
        if i == idx[0]:
            yield record
            idx = idx[1:]
            if not len(idx):
                break
        
def iterate_fasta(path, limit=None, random=False, random_state=None):
    if random_state is not None:
        np.random.seed(random_state)
    records = SeqIO.parse(path, 'fasta')
    if random:
        if limit is None:
            raise ValueError('Random can only be used together with limit')
        size = sum(1 for _ in SeqIO.parse(path, 'fasta'))
        limit = min(size, limit)
        print(f'Reading {limit:,} random sequences from: {path}')
        random_idx = np.array(sorted(np.random.choice(size, limit, replace=False)))
        records = iterate_fasta_index(records, random_idx)
    elif limit:
        records = itertools.islice(records, limit)

    return records

def save_input_fasta(inputs, labels, fasta_path, random=True, limit=None, random_state=None):
    i = 0
    index = []
    with open(fasta_path, 'w') as f:
        for path, label in zip(inputs, labels):
            for record in iterate_fasta(path, limit=limit, random=random, random_state=random_state):
                index.append(OrderedDict(
                    index=i,
                    label=label,
                    seq_id=record.id
                ))
                record.id = str(i)
                record.description = ''
                record.name = record.id
                SeqIO.write(record, f, 'fasta')
                i += 1

    index = pd.DataFrame(index)
    return index


def run_usearch(fasta_path, distance_path, maxdist, termdist=1.0):
    if not shutil.which('usearch'):
        raise UsumError('Missing usearch dependency on PATH. \nInstall it from: https://drive5.com/usearch/download.html')  
        
    cmd = [
        'usearch',
        '-calc_distmx', fasta_path,
        '-tabbedout', distance_path,
        '-maxdist', str(maxdist),
        '-termdist', str(termdist)
    ]
    print(f'Running USEARCH command: {" ".join(cmd)}')
    subprocess.check_output(cmd)
    
    
def load_sparse_dist_matrix(distance_path):
    dist_matrix = pd.read_csv(distance_path, header=None, sep='\t')
    
    diagonal = dist_matrix[0] == dist_matrix[1]
    row = np.concatenate([dist_matrix[0], dist_matrix[1][~diagonal]])
    col = np.concatenate([dist_matrix[1], dist_matrix[0][~diagonal]])
    data = 1 - np.concatenate([dist_matrix[2], dist_matrix[2][~diagonal]])

    dist_matrix = sparse.csr_matrix((data, (row, col)), dtype=np.float32)
    return 1 - dist_matrix.toarray()
    
    
def fit_umap(dist_matrix, random_state=None, neighbors=15, min_dist=0.1, spread=1.0):
    print(dist_matrix)
    reducer = umap.UMAP(
        n_neighbors=neighbors,
        random_state=random_state,
        min_dist=min_dist,
        spread=spread,
        metric='precomputed'
    )
    embedding = reducer.fit_transform(dist_matrix)
    return reducer, embedding

    
def fit_tsne(dist_matrix, random_state=None):
    reducer = TSNE(
        random_state=random_state,
        metric='precomputed'
    )
    embedding = reducer.fit_transform(dist_matrix)
    return reducer, embedding

# Taken from from umap.plot, adjusted for t-SNE
def plot_points(
    reducer_object,
    labels=None,
    values=None,
    theme=None,
    cmap="Blues",
    color_key=None,
    color_key_cmap="Spectral",
    background="white",
    width=800,
    height=800,
    show_legend=True,
):

    if theme is not None:
        cmap = _themes[theme]["cmap"]
        color_key_cmap = _themes[theme]["color_key_cmap"]
        background = _themes[theme]["background"]

    if labels is not None and values is not None:
        raise ValueError(
            "Conflicting options; only one of labels or values should be set"
        )

    points = reducer_object.embedding_

    if points.shape[1] != 2:
        raise ValueError("Plotting is currently only implemented for 2D embeddings")

    font_color = _select_font_color(background)

    dpi = plt.rcParams["figure.dpi"]
    fig = plt.figure(figsize=(width / dpi, height / dpi))
    ax = fig.add_subplot(111)

    if points.shape[0] <= width * height // 10:
        ax = _matplotlib_points(
            points,
            ax,
            labels,
            values,
            cmap,
            color_key,
            color_key_cmap,
            background,
            width,
            height,
            show_legend,
        )
    else:
        ax = _datashade_points(
            points,
            ax,
            labels,
            values,
            cmap,
            color_key,
            color_key_cmap,
            background,
            width,
            height,
            show_legend,
        )

    ax.set(xticks=[], yticks=[])

    return ax

