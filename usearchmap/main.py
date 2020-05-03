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

def main(argv=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", nargs='+', help="Input FASTA sequence file paths.")
    parser.add_argument("--labels", nargs='*', required=False, help="Input file labels.")
    parser.add_argument("--output", required=True, help="Output directory path.")
    parser.add_argument("--maxdist", type=float, required=True, help="USEARCH: Maximum distance which should be written.")
    parser.add_argument("--termdist", type=float, default=1.0, help="USEARCH: Identity threshold for terminating the calculation. This should be set higher than maxdist.")
    parser.add_argument("--neighbors", type=int, default=15, help="UMAP: The size of local neighborhood.")
    parser.add_argument("--theme", default='viridis', help="UMAP: Plot color theme.")
    parser.add_argument("--width", type=int, default=800, help="UMAP: Plot width in pixels.")
    parser.add_argument("--height", type=int, default=800, help="UMAP: Plot height in pixels.")
    
    options = parser.parse_args()

    if options.labels:
        labels = options.labels
        if len(options.labels) != len(options.inputs):
            raise ValueError('Labels have to be provided for each input file (--label a b c)')
    else:
        labels = [os.path.splitext(os.path.basename(path))[0] for path in inputs]

    usearchmap(
        inputs=options.inputs, 
        labels=labels, 
        output=options.output, 
        maxdist=options.maxdist, 
        termdist=options.termdist,
        neighbors=options.neighbors,
        theme=options.theme,
        width=options.width,
        height=options.height
    )

    
def usearchmap(inputs, labels, output, maxdist, termdist=1.0, neighbors=15, theme='viridis', width=800, height=800):
    print(f'> Creating directory: {output}')
    os.makedirs(output, exist_ok=True)
    fasta_path = os.path.join(output, 'input.fa')
    index_path = os.path.join(output, 'index.tsv')
    distance_path = os.path.join(output, 'distance.txt')
    png_path = os.path.join(output, 'umap.png')
    html_path = os.path.join(output, 'umap.html')
    
    print(f'> Writing FASTA records from {len(inputs)} paths...')
    index = save_input_fasta(inputs, labels, fasta_path)
    print(f'Saved {len(index):,} records to: {fasta_path}')
    
    index.to_csv(index_path, sep='\t', index=False)
    print(f'Saved index TSV to: {index_path}')
    
    print(f'> Creating sparse {len(index):,} x {len(index):,} distance matrix...')
    run_usearch(fasta_path, distance_path, maxdist=maxdist, termdist=termdist)
    dist_matrix = load_sparse_dist_matrix(distance_path)
    # with {dist_matrix.getnnx():,} non-null values
    print(f'Loaded {len(index):,} x {len(index):,} distance matrix ({dist_matrix.nbytes / 1024 / 1024} MB)')
    
    print(f'> Creating UMAP embedding with {neighbors} neighbors...')
    reducer = fit_umap(dist_matrix, neighbors=neighbors)
    
    print('> Drawing UMAP PNG...')
    ax = umap.plot.points(reducer, labels=index['label'], theme=theme, width=width, height=height);
    ax.figure.savefig(png_path)
    print(f'Saved UMAP PNG to: {png_path}')
    
    print('> Drawing interactive UMAP...')
    p = umap.plot.interactive(reducer, theme=theme, hover_data=index);
    bokeh.plotting.output_file(html_path)
    bokeh.plotting.save(p)
    print(f'Saved UMAP HTML to: {html_path}')
    
    print(f'Done. Saved to: {output}')

def save_input_fasta(inputs, labels, fasta_path):
    i = 0
    index = []
    with open(fasta_path, 'w') as f:
        for path, label in zip(inputs, labels):
            for r, record in enumerate(SeqIO.parse(path, 'fasta')):
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                # TODO
                if r % 100 == 0:
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
        print('Missing usearch dependency on PATH', file=sys.stderr)
        print('Install it from: https://drive5.com/usearch/download.html', file=sys.stderr)
        sys.exit(2)
        
    cmd = [
        'usearch',
        '-calc_distmx', fasta_path,
        '-tabbedout', distance_path,
        '-maxdist', str(maxdist),
        '-termdist', str(termdist)
    ]
    print(f'> Running USEARCH command: {" ".join(cmd)}')
    subprocess.check_output(cmd)
    
    
def load_sparse_dist_matrix(distance_path):
    dist_matrix = pd.read_csv(distance_path, header=None, sep='\t')

    row = dist_matrix[0]
    col = dist_matrix[1]
    data = 1-dist_matrix[2]

    dist_matrix = sparse.csr_matrix((data, (row, col)), dtype=np.float32)
    return (1 - dist_matrix.toarray())
    
    
def fit_umap(dist_matrix, neighbors=15):
    reducer = umap.UMAP(
        n_neighbors=neighbors,
        #metric='precomputed'
    )
    reducer.fit(dist_matrix)
    return reducer