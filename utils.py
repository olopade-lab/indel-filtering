import gzip
import io
import pandas as pd

def read_vcf(
        path,
        dtype={
            '#CHROM': str,
            'POS': int,
            'ID': str,
            'REF': str,
            'ALT': str,
            'QUAL': str,
            'FILTER': str,
            'INFO': str,
            'FORMAT': str,
            'NORMAL': str,
            'TUMOUR': str
        }
    ):
    if path.endswith('.gz'):
        with gzip.open(path, 'r') as f:
            lines = [l for l in f if not l.startswith(b'##')]
        return pd.read_csv(
            io.BytesIO(b''.join(lines)),
            dtype=dtype,
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
    else:
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
            ''.join(lines),
            dtype=dtype,
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
