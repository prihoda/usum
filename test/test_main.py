from usum import usum
import os

def get_test_file(path):
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(here, 'data', path)

def run_usearch(fasta_path, distance_path, **kwargs):
    vals = {}

    # add edges between even-even and odd-odd numbers
    for i in range(0, 2000, 2):
        for j in range(0, 2000, 2):
            if i < j:
                vals[(i, j)] = 0.3
                if i < 999 and j < 999:
                    vals[(i+1, j+1)] = 0.2

    # add edges between five-divisible numbers
    for i in range(0, 2000, 5):
        for j in range(0, 2000, 5):
            if i < j:
                vals[(i, j)] = 0.1

    for i in range(2000):
        vals[(i, i)] = 0

    with open(distance_path, 'w') as f:
        for (i, j), v in vals.items():
            f.write(f'{i}\t{j}\t{v}')
            f.write('\n')

def test_usum_umap(tmpdir, mocker):
    tmpdir = str(tmpdir)

    mock_usearch = mocker.patch('usum.main.run_usearch')
    mock_usearch.side_effect = run_usearch

    usum(
        inputs=[get_test_file('first.fa'), get_test_file('second.fa')], 
        output=tmpdir, 
        maxdist=0.3, 
        termdist=0.5
    )

    mock_usearch.assert_called_with(
        os.path.join(tmpdir, 'input.fa'), 
        os.path.join(tmpdir, 'distance.txt'), 
        maxdist=0.3, 
        termdist=0.5,
        threads=None,
        dbstep=None
    )

    assert os.path.exists(tmpdir)
    assert os.path.exists(os.path.join(tmpdir, 'umap.png'))

    
def test_usum_tsne(tmpdir, mocker):
    tmpdir = str(tmpdir)

    mock_usearch = mocker.patch('usum.main.run_usearch')
    mock_usearch.side_effect = run_usearch

    usum(
        inputs=[get_test_file('first.fa'), get_test_file('second.fa')], 
        output=tmpdir, 
        maxdist=0.3, 
        termdist=0.5,
        method='tsne'
    )

    mock_usearch.assert_called_with(
        os.path.join(tmpdir, 'input.fa'), 
        os.path.join(tmpdir, 'distance.txt'), 
        maxdist=0.3, 
        termdist=0.5,
        threads=None,
        dbstep=None
    )

    assert os.path.exists(tmpdir)
    assert os.path.exists(os.path.join(tmpdir, 'tsne.png'))
    