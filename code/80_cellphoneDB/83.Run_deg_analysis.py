 import pandas as pd
import sys
import os
print(sys.version)
out_path = '/geschwindlabshares/RexachGroup/Xia_Data/CCI/CellPhoneDB/output'

dx = "AD"

celltype = "majorC"

cpdb_file_path='/geschwindlabshares/RexachGroup/Xia_Data/CCI/CellPhoneDB/db/v4.1.0/cellphonedb.zip'
meta_file = 'meta.majorC.' + dx + ".meta_scATAC_clean.txt"
h5ad_file = "scATAC_clean." + dx + ".h5ad"
deg_file = "DEG." + dx + ".downsample.txt"
meta_file_path = os.path.join(out_path, meta_file)
counts_file_path=os.path.join(out_path, h5ad_file)
degs_file_path=os.path.join(out_path, deg_file)

metadata = pd.read_csv(meta_file_path, sep = '\t')
metadata.head(3)
metadata.shape

import anndata

adata = anndata.read_h5ad(counts_file_path)
adata.shape #cell x gene

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

flag = "DEG.majorC." + dx
flag

deg = pd.read_csv(degs_file_path,sep = '\t')
deg.clusters.unique()
metadata.majorC.unique()

from cellphonedb.src.core.methods import cpdb_degs_analysis_method

deconvoluted, means, relevant_interactions, significant_means = cpdb_degs_analysis_method.call(
  cpdb_file_path = cpdb_file_path,                            # mandatory: CellPhoneDB database zip file.
  meta_file_path = meta_file_path,                            # mandatory: tsv file defining barcodes to cell label.
  counts_file_path = counts_file_path,                        # mandatory: normalized count matrix.
  degs_file_path = degs_file_path,                            # mandatory: tsv file with DEG to account.
  counts_data = 'hgnc_symbol',                                # defines the gene annotation in counts matrix.
  threshold = 0.1,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
  result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
  separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
  debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
  output_path = out_path,                                     # Path to save results
  output_suffix = flag                                        # Replaces the timestamp in the output files by a user defined string in the  (default: None)
)



