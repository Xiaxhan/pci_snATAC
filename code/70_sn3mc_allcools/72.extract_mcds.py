import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

from ALLCools.mcds import MCDS
from ALLCools.clustering import tsne, significant_pc_test, filter_regions, lsi, binarize_matrix
from ALLCools.plot import *

import xarray as xr

mcds_path = 'DxPeak_mC.mcds/'
var_dim = 'bvFTD.odc.C7'
mcds = MCDS.open(mcds_path, obs_dim='cell', var_dim=var_dim)
mcds=MCDS.open(mcds_path)

mch_pattern='CNN'

mcad = mcds.get_score_adata(mc_type='CNN', quant_type='hypo-score')
mcad = mcds.get_score_adata(mc_type='CNN', quant_type='count')
df=mcad.T.to_df()
df.to_csv("test.csv")
