import pandas as pd
import sys
import os
print(sys.version)
out_path = 'CellPhoneDB/output'

dx = "bvFTD"

cpdb_file_path='CellPhoneDB/db/v4.1.0/cellphonedb.zip'
meta_file = 'meta.' + dx + ".meta_scATAC_clean.txt"
h5ad_file = "scATAC_clean." + dx + ".h5ad"
meta_file_path = os.path.join(out_path, meta_file)
counts_file_path=os.path.join(out_path, h5ad_file)

metadata = pd.read_csv(meta_file_path, sep = '\t')
metadata.head(3)
metadata.shape

import anndata

adata = anndata.read_h5ad(counts_file_path)
adata.shape #cell x gene

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

flag = dx

deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
  cpdb_file_path = cpdb_file_path,                 # mandatory: CellPhoneDB database zip file.
  meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
  counts_file_path = counts_file_path,             # mandatory: normalized count matrix.
  counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
  microenvs_file_path = None,       # optional (default: None): defines cells per microenvironment.
  iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
  threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
  threads = 4,                                     # number of threads to use in the analysis.
  debug_seed = 42,                                 # debug randome seed. To disable >=0.
  result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
  pvalue = 0.05,                                   # P-value threshold to employ for significance.
  subsampling = False,                             # To enable subsampling the data (geometri sketching).
  subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
  subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
  subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
  separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
  debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
  output_path = out_path,                          # Path to save results.
  output_suffix = flag                          # Replaces the timestamp in the output files by a user defined string in the  (default: None).
)
