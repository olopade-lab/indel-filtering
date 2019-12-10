from collections import OrderedDict
import os
import glob
import time

import pandas as pd
import numpy as np

import parsl
from parsl.app.app import python_app

from parsl.configs.htex_local import config
config.executors[0].max_workers = 14
# from gardner import config
parsl.load(config)

from utils import read_parquet

dtype = OrderedDict(
        [
            ('Chr', 'str'),
            ('Start', 'int64'),
            ('End', 'int64'),
            ('Ref', 'str'),
            ('Alt', 'str'),
            ('Func.refGene', 'str'),
            ('Gene.refGene', 'str'),
            ('GeneDetail.refGene', 'str'),
            ('ExonicFunc.refGene', 'str'),
            ('AAChange.refGene', 'str'),
            ('avsnp150', 'str'),
            ('ExAC_ALL', 'str'),
            ('ExAC_AFR', 'str'),
            ('ExAC_AMR', 'str'),
            ('ExAC_EAS', 'str'),
            ('ExAC_FIN', 'str'),
            ('ExAC_NFE', 'str'),
            ('ExAC_OTH', 'str'),
            ('ExAC_SAS', 'str'),
            ('AF', 'str'),
            ('AF_popmax', 'str'),
            ('AF_male', 'str'),
            ('AF_female', 'str'),
            ('AF_raw', 'str'),
            ('AF_afr', 'str'),
            ('AF_sas', 'str'),
            ('AF_amr', 'str'),
            ('AF_eas', 'str'),
            ('AF_nfe', 'str'),
            ('AF_fin', 'str'),
            ('AF_asj', 'str'),
            ('AF_oth', 'str'),
            ('non_topmed_AF_popmax', 'str'),
            ('non_neuro_AF_popmax', 'str'),
            ('non_cancer_AF_popmax', 'str'),
            ('controls_AF_popmax', 'str'),
            ('CLNALLELEID', 'str'),
            ('CLNDN', 'str'),
            ('CLNDISDB', 'str'),
            ('CLNREVSTAT', 'str'),
            ('CLNSIG', 'str'),
            ('Otherinfo', 'str'),
            ('unknown1', 'int64'),
            ('unknown2', 'int64'),
            ('CHROM', 'str'),
            ('POS', 'int64'),
            ('ID', 'str'),
            ('REF', 'str'),
            ('ALT', 'str'),
            ('QUAL', 'str'),
            ('FILTER', 'str'),
            ('INFO', 'str'),
            ('FORMAT', 'str'),
            ('NORMAL', 'str'),
            ('TUMOUR', 'str')
        ]
    )

@python_app
def get_annotated_df(paths, dtype):
    import os
    import time

    import pandas as pd

    samples = [os.path.basename(path).split('.')[0] for path in paths]
    result = None
    for sample, path in zip(samples, paths):
        start = time.time()
        df = pd.read_table(path, skiprows=1, names=dtype.keys())
        print('read {} in {:.1f} seconds'.format(path, time.time() - start))
        df = df[df['FILTER'] == 'PASS']
        df['SAMPLE'] = sample
        df = df.astype(dtype)

        if result is None:
            result = df
        else:
            result = pd.concat([result, df])

    output = '/scratch/t.cri.awoodard/indel-filtering/annotations_{}.pkl'.format('_'.join(samples))
    pd.to_pickle(result, output)

    return output

files = '/scratch/t.cri.awoodard/indel-filtering/annotations/*/*txt'
paths = glob.glob(files)
chunksize = 3
result = None

chunks = [paths[i:i + chunksize] for i in range(0, len(paths), chunksize)]
futures = [get_annotated_df(chunk, dtype) for chunk in chunks]

print('finished chunks')
dfs = [pd.read_pickle(f.result()) for f in futures]
df = pd.concat(dfs)

pon = read_parquet('/scratch/t.cri.awoodard/indel-filtering/pseudo_pon.parquet')
keys = list(pon.columns.values)
pon_index = pon.set_index(keys).index

dbsnp = read_parquet('/scratch/t.cri.awoodard/indel-filtering/All_20180423.parquet')
dbsnp_index = dbsnp.set_index(keys).index

genome1k = read_parquet('/scratch/t.cri.awoodard/indel-filtering/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.parquet')
genome1k_index = genome1k.set_index(keys).index

cosmic = read_parquet('/scratch/t.cri.awoodard/indel-filtering/CosmicCodingMuts.parquet')
cosmic_index = cosmic.set_index(keys).index

sample_index = df.set_index(keys).index

df['pseudo PON'] = 'pass'
df.loc[sample_index.isin(pon_index), 'pseudo PON'] = 'fail'

df['dbsnp'] = 'pass'
df.loc[sample_index.isin(dbsnp_index), 'dbsnp'] = 'fail'

df['genomes1k'] = 'pass'
df.loc[sample_index.isin(genome1k_index), 'genomes1k'] = 'fail'

df['cosmic'] = 'not included'
df.loc[sample_index.isin(cosmic_index), 'cosmic'] = 'included'

df.to_parquet('/scratch/t.cri.awoodard/indel-filtering/annotations.parquet')
