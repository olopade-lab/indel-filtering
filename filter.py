import os
import glob
import time

import pandas as pd
import numpy as np

import parsl
from parsl.app.app import python_app

from gardner import config

parsl.load(config)

files = '/scratch/t.cri.awoodard/indel-filtering/cgpPindel/*/*.vcf.gz'
chunksize = 2

@python_app
def get_filtered_df(paths):
    import os

    import pandas as pd
    import numpy as np

    from utils import read_vcf

    pon = pd.read_pickle('/home/t.cri.awoodard/indel-filtering/pseudo_pon.pkl')
    keys = list(pon.columns.values)
    pon_index = pon.set_index(keys).index

    race_map = pd.read_csv('/home/t.cri.awoodard/indel-filtering/race.csv')
    tcga_clinical = pd.read_csv('/home/t.cri.awoodard/indel-filtering/nationwidechildrens.org_clinical_patient_brca.txt', sep='\t')
    nigerian_clinical = pd.read_csv('/home/t.cri.awoodard/indel-filtering/WABCS_final_IHC_2019-06-14.csv')

    result = None
    columns = ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'SAMPLE', 'QUAL', 'RACE', 'ER status', 'PR status', 'HER2 status']
    samples = [path.split('/')[5] for path in paths]

    for path, sample in zip(paths, samples):
        print('starting {}'.format(path))
        df = read_vcf(path)
        df['SAMPLE'] = sample
        sample_index = df.set_index(keys).index

        df.loc[(sample_index.isin(pon_index)) & (df.FILTER == 'PASS'), 'FILTER'] = 'pseudo PON'

        df['RACE'] = 'nigerian'
        for row in race_map.iterrows(): # TODO: very inefficient-- eliminate this loop
            sample, race = row[1]
            df.loc[df.SAMPLE == sample, 'RACE'] = 'TCGA ' + race

        df['ER status'] = np.NaN
        df['PR status'] = np.NaN
        df['HER2 status'] = np.NaN
        for row in tcga_clinical.iterrows(): # TODO: very inefficient-- eliminate this loop
            df.loc[df.SAMPLE == row[1]['bcr_patient_barcode'], 'ER status'] = row[1]['er_status_by_ihc']
            df.loc[df.SAMPLE == row[1]['bcr_patient_barcode'], 'PR status'] = row[1]['pr_status_by_ihc']
            df.loc[df.SAMPLE == row[1]['bcr_patient_barcode'], 'HER2 status'] = row[1]['her2_status_by_ihc']

        for row in nigerian_clinical.iterrows(): # TODO: very inefficient-- eliminate this loop
            df.loc[df.SAMPLE == row[1]['NBCS Barcode'], 'ER status'] = row[1]['ER Status']
            df.loc[df.SAMPLE == row[1]['NBCS Barcode'], 'PR status'] = row[1]['PR Status']
            df.loc[df.SAMPLE == row[1]['NBCS Barcode'], 'HER2 status'] = row[1]['HER2 Status by IHC']

        if result is None:
            result = df[columns]
        else:
            result = pd.concat([result, df[columns]])
    output = '/scratch/t.cri.awoodard/indel-filtering/filtered_df_{}.pkl'.format('_'.join(samples))
    pd.to_pickle(result, output)

    for path in paths:
        os.unlink(path)

    return output

paths = glob.glob(files)
chunks = [paths[i:i + chunksize] for i in range(0, len(paths), chunksize)]

futures = [get_filtered_df(chunk) for chunk in chunks]

dfs = [pd.read_pickle(f.result()) for f in futures]
df = pd.concat(dfs)

df = df.sample(frac=1).reset_index(drop=True)
print('finished shuffling df')

df.to_parquet('/scratch/t.cri.awoodard/indel-filtering/cgp_pindel_filtered_shuffled_hormone_status.parquet')
