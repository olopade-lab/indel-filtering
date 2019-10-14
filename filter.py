import os
import glob

import pandas as pd

from multiprocessing import Pool

files = '/scratch/t.cri.awoodard/indel-filtering/cgpPindel/*/*.vcf.gz'
num_workers = 20
chunksize = 5

def get_filtered_df(paths):
    print('starting')
    from utils import read_vcf
    import pandas as pd

    pon = pd.read_pickle('pseudo_pon.pkl')
    keys = list(pon.columns.values)
    pon_index = pon.set_index(keys).index

    race_data = pd.read_csv('/home/t.cri.awoodard/indel-filtering/race.csv')

    result = None
    columns = ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'SAMPLE', 'QUAL', 'RACE']
    samples = [path.split('/')[5] for path in paths]

    for path, sample in zip(paths, samples):
        print('starting {}'.format(path))
        df = read_vcf(path)
        df['SAMPLE'] = sample
        sample_index = df.set_index(keys).index

        df.loc[(sample_index.isin(pon_index)) & (df.FILTER == 'PASS'), 'FILTER'] = 'pseudo PON'

        df['RACE'] = 'nigerian'
        for row in race_data.iterrows():
            sample, race = row[1]
            df.loc[df.SAMPLE == sample, 'RACE'] = 'TCGA ' + race

        if result is None:
            result = df[columns]
        else:
            result = pd.concat([result, df[columns]])
    output = '/scratch/t.cri.awoodard/indel-filtering/filtered_df_{}.pkl'.format('_'.join(samples))
    pd.to_pickle(result, output)

paths = glob.glob(files)
chunks = [paths[i:i + chunksize] for i in range(0, len(paths), chunksize)]

pool = Pool(processes=num_workers)
print('started pool')

pool.map(get_filtered_df, chunks)
print('finished filtering dfs')

dfs = [pd.read_pickle(path) for path in glob.glob('/scratch/t.cri.awoodard/indel-filtering/filtered_df_*.pkl')]

df = pd.concat(dfs)
print('finished concatenating dfs')

df = df.sample(frac=1).reset_index(drop=True)
print('finished shuffling df')

df.to_parquet('/scratch/t.cri.awoodard/indel-filtering/cgp_pindel_filtered_shuffled.parquet')
