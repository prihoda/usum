from setuptools import setup

setup(
    name='usearchmap',
    packages=['usearchmap'],
    install_requires=[
        'biopython',
        'umap-learn[plot]'
    ],
    entry_points={
        'console_scripts': ['usearchmap = usearchmap.main:main']
    }
)

