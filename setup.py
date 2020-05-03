from setuptools import setup
import os

about = {}
# Read version number from deepbgc.__version__.py (see PEP 396)
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'usum', '__version__.py'), encoding='utf-8') as f:
    exec(f.read(), about)

# Read contents of readme file into string
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='usum',
    packages=['usum'],
    version=about['__version__'],
    description='USUM: Plotting sequence similarity using USEARCH & UMAP',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='David Příhoda',
    author_email='david.prihoda@gmail.com',
    license='MIT',
    python_requires=">=3.6",
    keywords='dna, protein, sequence, similarity, umap, usearch, uclust, plot',
    url='https://github.com/prihoda/usum',
    install_requires=[
        'biopython',
        'umap-learn[plot]'
    ],
    entry_points={
        'console_scripts': ['usum = usum.main:main']
    }
)

