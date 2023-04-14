import os
import pathlib as pl
base = pl.Path(os.getcwd())

mklist = []
mklist.append(base / 'logs')
mklist.append(base / 'raw_data')
mklist.append(base / 'raw_data' / 'aux_scenic')

for stage in ['seurat', 'scenic']:
    mklist.append(base / 'output_data' / stage)
    mklist.append(base / 'figures' / stage)
for each in mklist:
    # print(each)
    each.mkdir(parents=True, exist_ok=True)