import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from ALLCools.mcds import MCDS
from ALLCools.clustering import tsne, significant_pc_test, filter_regions, lsi, binarize_matrix
from ALLCools.plot import *
import xarray as xr
import re
from os.path import exists
out_path = 'hyperScore'

mcds_path = 'dar_mC_subC'

dxlist = open("list", "r") # dar files with >1 dar in thatfile


for item in dxlist:
    group = re.findall(r"(.+?)\.sort", item)[0]
    groupdir = mcds_path + group
    if exists(groupdir):
        print(group)
        print(groupdir)
        mcds = MCDS.open(groupdir, obs_dim='cell', var_dim=group)

        mcad = mcds.get_score_adata(mc_type='CNN', quant_type='hyper-score')
        df=mcad.T.to_df()
        outfile="hyperScore." + group + ".csv"
        outfile = out_path + "/" + outfile
        df.to_csv(outfile)

        print(out_path + "/" + outfile)
    else:
        print(group)
