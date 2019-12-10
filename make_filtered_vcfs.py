import os
import glob

import pandas as pd
import numpy as np

import parsl
from parsl.app.app import python_app

from gardner import config

parsl.load(config)

base = '/scratch/t.cri.awoodard/indel-filtering/filtered_cgpPindel'
files = '/scratch/t.cri.awoodard/indel-filtering/cgpPindel/*/*.vcf.gz'
chunksize = 20

@python_app(cache=True)
def make_filtered_vcfs(paths):
    import os
    import time
    import gzip
    import io
    from pathlib import Path

    import pandas as pd
    import numpy as np

    from utils import read_parquet

    keys = ['CHROM', 'POS', 'REF', 'ALT']

    start = time.time()
    dbsnp = read_parquet('/scratch/t.cri.awoodard/indel-filtering/All_20180423.parquet')
    dbsnp_index = dbsnp.set_index(keys).index

    genome1k = read_parquet('/scratch/t.cri.awoodard/indel-filtering/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.parquet')
    genome1k_index = genome1k.set_index(keys).index

    cosmic = read_parquet('/scratch/t.cri.awoodard/indel-filtering/CosmicCodingMuts.parquet')
    cosmic_index = cosmic.set_index(keys).index
    print('finished indexing in {:.0} seconds'.format(time.time() - start))

    for path in paths:
        print('starting {}'.format(path))

        patient_id, filename = path.split('/')[-2:]
        os.makedirs(os.path.join(base, patient_id), exist_ok=True)
        if os.path.exists(os.path.join(base, patient_id, 'done')):
            continue

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

        with gzip.open(path, 'r') as f:
            lines = f.readlines()
            header = [l for l in lines if l.startswith(b'##')]
            body = [l for l in lines if not l.startswith(b'##')]
            chunks = pd.read_csv(
                io.BytesIO(b''.join(body)),
                dtype=dtype,
                sep='\t',
                iterator=True,
                chunksize=1000
            )
            # iterate through data, discarding non-pass as we go, to save memory
            df = pd.concat([chunk[chunk['FILTER'] == 'PASS'] for chunk in chunks]).rename(columns={'#CHROM': 'CHROM'})

        sample_index = df.set_index(keys).index

        start = time.time()
        dbsnp_pass = (~sample_index.isin(dbsnp_index) | sample_index.isin(cosmic_index))
        genome1k_pass = (~sample_index.isin(genome1k_index) | sample_index.isin(cosmic_index))
        print('finished filtering in {:.0f} seconds'.format(time.time() - start))

        with gzip.open(os.path.join(base, patient_id, 'pass_dbSNP_1kgenomes_' + filename), 'wb') as f:
            f.write(b''.join(header))
            f.write(b'#' + df[dbsnp_pass & genome1k_pass].to_csv(index=False, sep='\t').encode())
        with gzip.open(os.path.join(base, patient_id, 'fail_dbSNP_1kgenomes_' + filename), 'wb') as f:
            f.write(b''.join(header))
            f.write(b'#' + df[~dbsnp_pass | ~genome1k_pass].to_csv(index=False, sep='\t').encode())

        Path(os.path.join(base, patient_id, 'done')).touch()

paths = glob.glob(files)
chunks = [paths[i:i + chunksize] for i in range(0, len(paths), chunksize)]
futures = [make_filtered_vcfs(chunk) for chunk in chunks]
parsl.wait_for_current_tasks()

print('finished!')

